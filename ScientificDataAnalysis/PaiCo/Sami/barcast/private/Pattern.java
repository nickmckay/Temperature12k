import java.util.*;

public class Pattern {

    public int item = -1;
    
    public HashMap<Integer, Pattern> children = new HashMap<Integer,Pattern>();
    public LinkedList<Integer> indices = new LinkedList<Integer>();
    
    public Pattern() {
    }
    
    public Pattern(int item) {
        this.item = item;
    }
    
    public int findCommon(int[] items1, int[] items2) {
        Pattern p1 = this;
        Pattern p2 = this;
        for (int i = 0; i < items1.length & p1 != null; i++) {
            p1 = p1.children.get(items1[i]);
        }
        for (int i = 0; i < items2.length & p2 != null; i++) {
            p2 = p2.children.get(items2[i]);
        }
        
        if (p1 != null & p2 != null) {
            Iterator i1 = p1.indices.iterator();
            Iterator i2 = p2.indices.iterator();
            int v1 = -2;
            int v2 = -1;
            while (v1 != v2) {
                if (v1 < v2) {
                    if (!i1.hasNext()) break;
                    v1 = ((Integer)i1.next()).intValue();
                } else if (v2 < v1) {
                    if (!i2.hasNext()) break;
                    v2 = ((Integer)i2.next()).intValue();
                } 
            }     
            if (v1 == v2)
                return v1;
        }
        return -1;
    }
    
    /**
     * Returns the stored patterns in the prefix tree and the mapping
     * from indices to patterns. The last element in the resulting array
     * is the mapping.
     */
    public int[][] getPatterns() {
        LinkedList<LinkedList<Integer>> patterns = new LinkedList<LinkedList<Integer>>();
        Stack<Pattern> stack = new Stack<Pattern>();
        for(Pattern pat: children.values()) {
            stack.push(pat);
        }
        LinkedList<Integer> items = new LinkedList<Integer>();
        HashMap<Integer, Integer> map = new HashMap<Integer,Integer>();
        
        // Depth-firsth search to find all patterns
        Integer maxIndex = -1;
        while (!stack.empty()) {
            Pattern top = stack.pop();
            if (top == null) {
                items.removeLast();
            } else {
                items.add(top.item);
                if (top.indices.size() > 0) {
                    patterns.add(new LinkedList<Integer>(items));
                    for(Integer i: top.indices) {
                        map.put(i, patterns.size());
                        if (maxIndex < i)
                            maxIndex = i;
                    }
                }
                
                // Null signifies depth step
                stack.push(null);
                for (Pattern pat: top.children.values()) {
                    stack.push(pat);
                }
            }
        }
        
        int result[][] = new int[patterns.size()+1][];
        int index = 0;
        // Store found patterns to result
        for(LinkedList<Integer> pat: patterns) {
            result[index] = new int[pat.size()];
            int j = 0;
            for(Integer i: pat) {
                result[index][j++] = i.intValue();
            }
            index++;
        }
        
        // Store also the mapping from indices to patterns
        result[index] = new int[maxIndex];
        for (Integer key: map.keySet()) {
            result[index][key-1] = map.get(key);
        }
        
        return result;        
    }
    
    /**
     * Stores index to the prefix tree node of items, extending the tree 
     * necessary
     * @param items pattern to store
     * @param index index to store
     */
    public void store(int[] items, int index) {
        Pattern p = this;
        for (int i = 0; i < items.length; i++) {
            Pattern np = p.children.get(items[i]);
            if (np == null) {                
                np = new Pattern(items[i]);                
                p.children.put(items[i],np);
            }
            p = np;
        }
        p.indices.add(index);
    }

    public String toString(String prefix) {
        String s= "";
        Iterator i = children.keySet().iterator();
        while (i.hasNext()) {
            Integer key = (Integer)i.next();
            s += children.get(key).toString(prefix + key + " ");
        }
        
        if (indices.size() > 0) {
            Iterator j = indices.iterator();
            String t = "";
            while (j.hasNext()) {
                Integer k = (Integer)j.next();
                if (j.hasNext()) 
                    t += k + " ";
                else 
                    t += k;
            }
            s += "(" + prefix + ")=[" + t +"]\n"; 
        }
        return s;
    }    
    
    public String toString() {
        return toString("");
    }    
};