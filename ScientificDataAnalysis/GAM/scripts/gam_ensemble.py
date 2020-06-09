import pyleogrid as pg
import numpy as np
import pygam
import datetime as dt
import psyplot.data as psyd
import pandas as pd
import xarray as xr
from sklearn.neighbors import BallTree
import distributed


def is_between(da, left, right):
    return (da >= left) & (da <= right)


def get_ensemble_overlap(ens1, ens2):
    return (is_between(ens1, ens2.min(), ens2.max()),
            is_between(ens2, ens1.min(), ens1.max()))


class GAMEnsemble(pg.Ensemble):
    """A :class:`pyleogrid.Ensemble` using pseudo-gridding with GAM"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # use psyplot base to keep coords of unstructured datasets
        if self.target is not None:
            self.output_ds = self.target.psy.base

    def find_gridcells(self, target):
        conf = self.config

        full_coords = self.input_data.reset_coords()[
            [conf.lon, conf.lat]].to_dataframe()
        coords = full_coords.drop_duplicates().copy()

        target_lon = target.psy.decoder.get_x(target.psy[0])
        target_lat = target.psy.decoder.get_y(target.psy[0])

        # load grid cell centers
        tree = BallTree(np.r_[target_lon.values[None],
                              target_lat.values[None]].T)
        indices = tree.query(coords, return_distance=False)[:, 0]

        lonname = target_lon.name
        if lonname == conf.lon:
            lonname = conf.lon + "_grid"
        latname = target_lat.name
        if latname == conf.lat:
            latname = conf.lat + "_grid"

        coords[lonname] = target_lon.values[indices]
        coords[latname] = target_lat.values[indices]

        full_coords = full_coords.merge(coords, on=[conf.lon, conf.lat],
                                        how='left')
        agedim = self.input_data[conf.age].dims[0]
        self.input_data[lonname] = (agedim, full_coords[lonname])
        self.input_data[latname] = (agedim, full_coords[latname])

        self.input_data = self.input_data.set_coords(
            [lonname, latname]).sortby(
                [latname, lonname, conf.ds_id, conf.age])

        return self.input_data[lonname], self.input_data[latname]

    def predict(
            self, age=None, climate=None, target=None,
            size=None, client=None, quantiles=None, return_gam=False,
            return_time=False, return_counts=False,
            **kwargs):

        conf = self.config

        if target is None:
            target = self.target

        if target is None:
            raise ValueError("No target elevation data supplied!")

        clon, clat = self.find_gridcells(target)

        if age is None:
            age = self.input_data.get(conf.age + '_ensemble',
                                      self.input_data[conf.age])
        if climate is None:
            climate = self.input_data.get(conf.climate + '_ensemble',
                                          self.input_data.get(conf.climate))
            if climate is None:
                raise ValueError("No climate data found!")

        data = self.input_data

        agedim = data[conf.age].dims[0]

        data[conf.climate + '_ensemble'] = climate
        data[conf.age + '_ensemble'] = age

        output = target.psy.base.copy()

        time = target.dims[0]  # time dimension
        space = target.dims[1]  # spatial dimension
        nt = output.dims[time]
        ns = output.dims[space]

        encoding = dict(coordinates=f'{clon.name} {clat.name}')

        output[conf.climate] = target.copy(
            data=np.full(target.shape, np.nan))
        if quantiles is not None:
            output[conf.climate + '_quantiles'] = xr.Variable(
                (time, space, 'quantile'),
                np.full((nt, ns, len(quantiles)), np.nan),
                encoding=encoding)
            output['quantile'] = ('quantile', quantiles)
        if size is not None:
            output[conf.climate + '_samples'] = xr.Variable(
                (time, conf.ensemble, space), np.full((nt, size, ns), np.nan),
                encoding=encoding)
        if return_gam:
            output[conf.climate + '_model'] = xr.Variable(
                space, np.full(ns, object, dtype=object),
                encoding=encoding)
        if return_time:
            output['time_needed'] = xr.Variable(
                space, np.zeros(ns, int), encoding=encoding)
        if return_counts:
            output['nsamples'] = target.copy(
                data=np.full(target.shape, np.nan))

        time = target.coords[target.dims[0]].values

        kwargs.update(dict(time=time, quantiles=quantiles, size=size,
                           return_gam=return_gam, return_counts=return_counts))

        target_idx = pd.MultiIndex.from_arrays(
            [target.psy.decoder.get_x(target.psy[0]).values,
             target.psy.decoder.get_y(target.psy[0]).values])

        kwargs['conf'] = conf

        coords = pd.MultiIndex.from_arrays(
            [clon.values, clat.values])

        def iterate():
            locs = map(coords.get_loc, coords.drop_duplicates())
            if client is not None:
                futures = client.map(
                    self._predict_cell,
                    [(sl, data.isel(**{agedim: sl})) for sl in locs],
                    **kwargs)
                for future, result in distributed.as_completed(
                        futures, with_results=True):
                    futures.remove(future)
                    yield result
            else:
                for sl in locs:
                    yield self._predict_cell((sl, data.isel(**{agedim: sl})),
                                             **kwargs)

        for ret in pg.utils.log_progress(
                iterate(), len(coords.drop_duplicates())):
            if ret is None:
                continue
            key = ret[0]
            i = target_idx.get_loc(coords[key][0])
            try:
                i = i.start
            except AttributeError:
                pass
            output[conf.climate][:, i] = ret[1]
            j = 2
            if quantiles is not None:
                output[conf.climate + '_quantiles'][:, i, :] = ret[j]
                j += 1
            if size is not None:
                output[conf.climate + '_samples'][:, :, i] = ret[j]
                j += 1
            if return_gam:
                output[conf.climate + '_model'][i] = ret[j]
                j += 1
            if return_time:
                output['time_needed'][i] = ret[j]
                j += 1
            if return_counts:
                output['nsamples'][:, i] = ret[j]
                j += 1

        return output

    @staticmethod
    def _predict_cell(keyds, conf, time, quantiles=None, size=None,
                      align=True, anomalies=True, min_overlap=100,
                      return_gam=False, max_time_diff=200, return_time=False,
                      modern_young=-50, modern_old=-20, ensemble_cls=None,
                      use_gam=True, return_counts=False):

        t0 = dt.datetime.now()

        key, ds = keyds

        agedim = ds[conf.age].dims[0]

        if ensemble_cls is None:
            ensemble_cls = GAMEnsemble

        if align:
            ds = ensemble_cls._align_ensembles(
                ds, conf, min_overlap=min_overlap)
        else:
            ds['aligned'] = (agedim, np.zeros(ds.dims[agedim], bool))

        if anomalies:
            ds = ensemble_cls._compute_anomaly(
                ds, conf, modern_young, modern_old)

        if use_gam:
            ret = ensemble_cls._predict_gam(
                ds, conf, time, quantiles, size, return_gam,
                return_counts, max_time_diff)
        else:
            ret = ensemble_cls._binned_mean(
                ds, conf, time, quantiles, size, return_counts, max_time_diff)
        if ret is None:
            return

        if return_time:
            ret = ret + ((dt.datetime.now() - t0).total_seconds(), )

        return (key, ) + ret

    @staticmethod
    def _predict_gam(ds, conf, time, quantiles=None, size=None,
                     return_gam=False,  return_counts=False,
                     max_time_diff=200):
        # insert 0s for every timeseries in the ensemble for the reference
        # period at -35 BP (1985)

        climate = conf.climate + '_ensemble'
        age = conf.age + '_ensemble'

        x = ds[age].values.ravel()
        y = ds[climate].values.ravel()

        mask = (~np.isnan(x)) & (~np.isnan(y))
        if not mask.any():
            return
        else:
            x = x[mask]
            y = y[mask]

        gam = pygam.LinearGAM(pygam.s(0)).gridsearch(
            x[:, np.newaxis], y, progress=False)

        time = np.asarray(time)

        ret = (gam.predict(time), )

        if quantiles is not None:
            ret = ret + (gam.prediction_intervals(time, quantiles=quantiles), )
        if size is not None:
            ret = ret + (gam.sample(
                x[:, np.newaxis], y, sample_at_X=time, n_draws=size).T, )
        if return_counts:
            tree = BallTree(ds[age].values.ravel()[:, np.newaxis])
            counts = tree.query_radius(time[:, np.newaxis], return_counts,
                                       count_only=True)
            ret = ret + (counts, )

        # look how many samples in the ensemble fall into the `max_time_diff`
        # time interval around the predicted time
        tree = BallTree(ds[age].values.ravel()[:, np.newaxis])
        counts = tree.query_radius(time[:, np.newaxis], max_time_diff,
                                   count_only=True)

        idx = counts < 100
        if idx.any():
            for arr in ret:
                arr[idx] = np.nan

        if return_gam:
            return ret + (gam, )
        else:
            return ret

    def _binned_mean(ds, conf, time, quantiles=None, size=None,
                     return_counts=False, max_time_diff=200):
        climate = conf.climate + '_ensemble'
        age = conf.age + '_ensemble'
        agedim = ds[conf.age].dims[0]
        ens = conf.ensemble

        def bootstrap_mean(da):
            resampler = np.random.randint(0, len(da), (len(da), size))
            return xr.DataArray(da.values[resampler].mean(axis=0),
                                dims=(ens, ))

        ds = ds.stack(**{age + ens: (agedim, ens)})

        time = pd.Index(time)
        tree = BallTree(time[:, np.newaxis])
        ind, dists = tree.query_radius(
            ds[age].values[:, np.newaxis], max_time_diff,
            return_distance=True, sort_results=True)

        miss = len(time)

        ind = np.array([t[0] if t.size else miss for t in ind])

        grouper = ds[age + ens].copy(data=np.r_[time, [np.nan]][ind])

        grouped = ds[climate].groupby(grouper)

        mask = grouped.count() > 100

        ret = (grouped.mean().where(mask), )
        if quantiles is not None:
            ret = ret + (grouped.quantile(quantiles).where(mask), )
        if size is not None:
            ret = ret + (grouped.apply(bootstrap_mean).where(mask), )
        if return_counts:
            tree = BallTree(ds[age].values.ravel()[:, np.newaxis])
            counts = tree.query_radius(time[:, np.newaxis], return_counts,
                                       count_only=True)
            ret = ret + (xr.DataArray(counts, dims=ret[0].dims[0],
                                      coords={ret[0].dims[0]: time}), )

        return tuple(arr.reindex({age + ens: time}).values for arr in ret)

    def align_ensembles(self, target=None, client=None, min_overlap=100,
                        inplace=True):
        """Align the ensembles per grid cell in the given target"""

        conf = self.config

        if target is None:
            target = self.target

        if target is None:
            raise ValueError("No target elevation data supplied!")

        clon, clat = self.find_gridcells(target)

        data = self.input_data

        agedim = data[conf.age].dims[0]

        coords = pd.MultiIndex.from_arrays(
            [clon.values, clat.values])

        kwargs = dict(min_overlap=min_overlap, conf=conf)

        def iterate():
            locs = map(coords.get_loc, coords.drop_duplicates())
            if client is not None:
                futures = client.map(
                    self._align_ensembles,
                    [data.isel(**{agedim: sl}) for sl in locs],
                    **kwargs)
                for future, result in distributed.as_completed(
                        futures, with_results=True):
                    futures.remove(future)
                    yield result
            else:
                for sl in locs:
                    yield self._align_ensembles(data.isel(**{agedim: sl}),
                                                **kwargs)

        datasets = []
        for ret in pg.utils.log_progress(
                iterate(), len(coords.drop_duplicates())):
            datasets.append(ret)

        if len(datasets) > 1:
            ds = xr.concat(
                datasets, agedim, compat='override', coords='all')
        else:
            ds = datasets[0]

        if inplace:
            self.input_data = ds
        else:
            return ds

    @staticmethod
    def _align_ensembles(ds, conf, min_overlap=100):

        climate = conf.climate + '_ensemble'

        ds = ds.load()

        ds['aligned'] = xr.zeros_like(ds[conf.age], dtype=bool)
        ds['alignment_base'] = xr.zeros_like(ds[conf.age], dtype=bool)
        ds['unaligned_' + climate] = ds[climate].copy()
        groups = dict(ds.groupby(conf.ds_id))
        lengths = {ts_id: _ds[conf.age].max() -
                   _ds[conf.age].min()
                   for ts_id, _ds in groups.items()}
        max_id = max(lengths, key=lengths.get)
        aligned = groups.pop(max_id)
        aligned['alignment_base'][:] = True
        found_overlap = True

        agedim = ds.age.dims[0]

        ageens = conf.age + '_ensemble'

        while groups and found_overlap:
            found_overlap = False
            for ts_id in groups:
                m1, m2 = get_ensemble_overlap(aligned[ageens],
                                              groups[ts_id][ageens])
                if m1.sum() > min_overlap and m2.sum() > min_overlap:
                    found_overlap = True
                    ts_ds = groups.pop(ts_id)
                    diff = (
                        aligned[climate].values[m1.values].mean() -
                        ts_ds[climate].values[m2.values].mean())
                    ts_ds[climate] = ts_ds[climate] + diff
                    aligned = xr.concat([aligned, ts_ds], dim=agedim)
                    break

        aligned['aligned'][:] = True
        if groups:
            aligned = xr.concat([aligned] + list(groups.values()), dim=agedim)
        return aligned

    def compute_anomalies(
            self, target=None, client=None, modern_young=-50, modern_old=-20,
            inplace=True):
        """Align the ensembles per grid cell in the given target"""

        conf = self.config

        if target is None:
            target = self.target

        if target is None:
            raise ValueError("No target elevation data supplied!")

        clon, clat = self.find_gridcells(target)

        data = self.input_data

        agedim = data[conf.age].dims[0]

        coords = pd.MultiIndex.from_arrays(
            [clon.values, clat.values])

        kwargs = dict(modern_young=modern_young, modern_old=modern_old,
                      conf=conf)

        def iterate():
            locs = map(coords.get_loc, coords.drop_duplicates())
            if client is not None:
                futures = client.map(
                    self._compute_anomaly,
                    [data.isel(**{agedim: sl}) for sl in locs],
                    **kwargs)
                for future, result in distributed.as_completed(
                        futures, with_results=True):
                    futures.remove(future)
                    yield result
            else:
                for sl in locs:
                    yield self._compute_anomaly(data.isel(**{agedim: sl}),
                                                **kwargs)

        datasets = []
        for ret in pg.utils.log_progress(
                iterate(), len(coords.drop_duplicates())):
            datasets.append(ret)

        if len(datasets) > 1:
            ds = xr.concat(
                datasets, agedim, compat='override', coords='all')
        else:
            ds = datasets[0]

        if inplace:
            self.input_data = ds
        else:
            return ds

    @staticmethod
    def _compute_anomaly(ds, conf, modern_young=-50, modern_old=-20,
                         modern='modern'):

        agedim = ds[conf.age].dims[0]

        climate = conf.climate + '_ensemble'
        age = conf.age + '_ensemble'
        anom_ref = conf.climate + '_anomaly_ref'

        ds = ds.load()

        ds[anom_ref] = xr.zeros_like(ds[conf.age])

        # remove the anomaly for the aligned data
        if 'aligned' in ds:
            aligned_ids = np.where(ds.aligned)[0]
            aligned = ds.isel(age=aligned_ids)
            nonaligned = ds[conf.ds_id].where(~ds.aligned, drop=True)

            # define modern via worldclim reference: 1970-2000
            is_modern = is_between(
                aligned[age], modern_young, modern_old).values
            if is_modern.sum() > 100:
                mean = aligned[climate].values[is_modern].mean()
            elif aligned.datum[aligned.alignment_base][0] == 'anom':
                mean = 0
            elif modern in aligned:
                mean = aligned[modern][aligned.alignment_base][0].values
            else:
                mean = np.nan
            ds[climate].values[:, aligned_ids] = (
                ds[climate].values[:, aligned_ids] - mean)
            ds[anom_ref].values[aligned_ids] = mean

        else:
            nonaligned = ds[conf.ds_id]

        # remove anomaly for non-aligned data
        for ts_id in np.unique(nonaligned):
            ids = np.where(ds[conf.ds_id] == ts_id)[0]
            sub = ds.isel(**{agedim: ids})
            is_modern = is_between(
                sub[age], modern_young, modern_old).values
            if is_modern.sum() >= 100:
                mean = sub[climate].values[is_modern].mean()
            elif sub.datum[sub.alignment_base][0] == 'anom':
                mean = 0
            else:
                mean = sub[modern].values[np.newaxis]
            ds[climate].values[:, ids] = ds[climate].values[:, ids] - mean
            ds[anom_ref].values[ids] = mean

        return ds


def test_align_ensembles():

    longest = pd.DataFrame({
        'TSid': [1] * 6,
        'age': np.arange(6),
        'temperature': np.arange(6)})

    overlapping = pd.DataFrame({
        'TSid': [2] * 3,
        'age': np.arange(2, 5) + 0.5,
        'temperature': np.arange(2, 5) - 2})

    not_overlapping = pd.DataFrame({
        'TSid': [3] * 3,
        'age': [7, 8, 9],
        'temperature': [-1] * 3})

    combined = pd.concat([longest, overlapping, not_overlapping],
                         ignore_index=True)

    combined['lat'] = 0
    combined['lon'] = 0
    target = xr.DataArray(
        [0], coords={'clon': xr.Variable(('cell', ), [0]),
                     'clat': xr.Variable(('cell', ), [0])},
        name='temperature',
        dims=('cell', ), attrs={'coordinates': 'clon clat'})

    combined['age_unc'] = 0.1

    ensemble = GAMEnsemble(combined.set_index('age'), ds_id='TSid',
                           climate='temperature', target=target)
    age_ensemble = ensemble.sample_ages(size=1000)
    ensemble.input_data['temperature_ensemble'] = age_ensemble.copy(
        data=np.tile(combined.temperature.values[None], (1000, 1)))

    aligned = ensemble.align_ensembles(inplace=False)
    assert aligned.aligned.where(aligned.TSid == 1, drop=True).all()
    assert aligned.aligned.where(aligned.TSid == 2, drop=True).all()
    assert not aligned.aligned.where(aligned.TSid == 3, drop=True).any()
    shifted_ds = aligned.where(aligned.TSid == 2, drop=True)
    assert np.allclose(
        shifted_ds.temperature_ensemble.values,
        shifted_ds.unaligned_temperature_ensemble.values + 2.5,
        atol=0.2)
