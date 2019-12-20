// COPYRIGHT
//   2011, Sami Hanhijärvi <sami.hanhijarvi@iki.fi>
//
// LICENSE
//   This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
// 
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <map>
#include <stack>
#include <list>
#include <vector>
#include <sstream>
#include <string.h>

// #include "MersenneTwister.h"

// Fix redefinition conflict between Matlab R2010a and MS Visual Studio 2010
#ifdef _CHAR16T
#define CHAR16_T
#endif
#include "mex.h"

#define ABS(x) (((x)>0)?(x):(-x))
inline double SQR(double x) { return x*x;}

using namespace std;

class Node {
public:    
    
    int item;    
    map<unsigned short int, Node *> children;
    list<unsigned int> indices;
    
    Node(int newItem): item(newItem) {
    }
    
    int findCommon(double* items1, unsigned int len1, double *items2, unsigned int len2) {
        Node *p1 = this;
        Node *p2 = this;
        for (unsigned int i = 0; i < len1; i++) {
            map<unsigned short int, Node *>::iterator iter = p1->children.find((unsigned short int)items1[i]);
            if (iter != p1->children.end())
                p1 = iter->second;
            else {
                p1 = 0;
                break;
            }
        }
        for (unsigned int i = 0; i < len2; i++) {
            map<unsigned short int, Node *>::iterator iter = p2->children.find((unsigned short int)items2[i]);
            if (iter != p2->children.end())
                p2 = iter->second;
            else {
                p2 = 0;
                break;
            }
        }
        
        if (p1 && p2) {
            list<unsigned int>::iterator i1 = p1->indices.begin();
            list<unsigned int>::iterator i2 = p2->indices.begin();
            unsigned int v1 = *i1;
            unsigned int v2 = *i2;
            while (v1 != v2) {
                if (v1 < v2) {
                    i1++;
                    if (i1 == p1->indices.end()) break;
                    v1 = *i1;
                } else if (v2 < v1) {
                    i2++;
                    if (i2 == p2->indices.end()) break;
                    v2 = *i2;
                } 
            }     
            if (v1 == v2)
                return v1;
        }
        return -1;
    }
    
    /**
     * Stores index to the prefix tree node of items, extending the tree 
     * necessary
     * @param items pattern to store
     * @param index index to store
     */
    void store(double* items, unsigned int len, unsigned int index) {
        Node *p = this;
        for (unsigned int i = 0; i < len; i++) {
            map<unsigned short int, Node *>::iterator iter = p->children.find((unsigned short int)items[i]);
            if (iter != p->children.end())
                p = iter->second;
            else {
                Node *np = new Node((unsigned int)items[i]);
                p->children[(unsigned int)items[i]]= np;
                p = np;
            }            
        }
        p->indices.push_back(index);
    }

    void toString(stringstream & out, list<unsigned short int> prefix) {
        
        if (item >= 0)
            prefix.push_back(item);
        
        if (indices.size() > 0) {
            out << '(';
            for (list<unsigned short int>::iterator i = prefix.begin(); i != prefix.end(); i++) {
                out << *i << ' ';
            }
            out << ")=[";
            for (list<unsigned int>::iterator i = indices.begin(); i != indices.end(); i++) {
                out << *i << ' ';
            }
            out << ']' << endl;
        }
        
        
        map<unsigned short int, Node *>::iterator i;
        for (i = children.begin(); i != children.end(); i++) {
            i->second->toString(out, prefix);
        }                        
    }    
    
    string toString() {
        stringstream out;
        list<unsigned short int> prefix;
        this->toString(out, prefix);
        return out.str();
    }    
};

void convert(Node * tree, vector<vector<int> > &lowArray, vector<vector<int> > &highArray) {
    stack<Node *> todo;
    list<int> prefix;
    for (map<unsigned short int, Node*>::iterator iter = tree->children.begin(); iter != tree->children.end(); iter++) {
        todo.push(iter->second);
    }
    while (todo.size() > 0) {
        Node * node = (Node *)todo.top();
        
        if (node->item  < 0) {
            // Remove the last item from prefix
            prefix.pop_back(); 
            // Delete the node as we don't need it anymore
            delete node;
            // Reduce the stack.
            todo.pop();
        } else {
            prefix.push_back(node->item);           
            
            // Add the prefix pattern to the arrays to positions indicated by 
            // the nodes indices.
            list<unsigned int>::iterator iter;
            for (iter = node->indices.begin(); iter != node->indices.end(); iter++) {
                if (highArray[*iter].size() == 0) {

                    // highArray stores the time patterns that have a larger calendar time than
                    // the lowArray time patterns. This is assured by the fact that patterns
                    // are stored in ascending index order in the prefix tree and when the tree is
                    // depth-first-searched in inverted order (as done here), the patterns are 
                    // iterated in descending order. As time index patterns are disjoint, this 
                    // means that high calendar patterns (large indices) are necessarily visited first.
                    highArray[*iter].resize(prefix.size());
                    list<int>::iterator pi;
                    int pos = 0;
                    for (pi = prefix.begin(); pi != prefix.end(); pi++) {
                        highArray[*iter][pos++] = *pi-1;
                    }
                } else {
                    if (lowArray[*iter].size() > 0) {
                        mexErrMsgTxt("Index collision.");
                    }
                    lowArray[*iter].resize(prefix.size());
                    list<int>::iterator pi;
                    int pos = 0;
                    for (pi = prefix.begin(); pi != prefix.end(); pi++) {
                        lowArray[*iter][pos++] = *pi-1;
                    }
                }
            }

            // Change the node item to -1 to mark it as a stop in the stack
            // to indicate the branch has been thoroughly visited and the 
            // prefix needs to be shortened. The node is also deleted.
            node->item = -1;
            
            // Add all children to the stack and iterate
            for (map<unsigned short int, Node *>::iterator ci = node->children.begin(); ci != node->children.end(); ci++) {
                todo.push(ci->second);
            }            
        }        
    }
}

void printUsage() {
    mexPrintf("paico_mex version 1.0\n");
    mexPrintf("USAGE:\n");
    mexPrintf("--Prefix tree--\n");
    mexPrintf("       paico_mex()                  Clears the prefix three and creates a new one.\n");
    mexPrintf(" ind = paico_mex(list1, list2)      Checks if they have the same index in the tree.\n");
    mexPrintf("       paico_mex(list, index)       Stores the prefix list with the given index.\n");
    mexPrintf(" str = paico_mex()                  Returns the prefix tree as a string for printing purposes.\n");
    mexPrintf("--Index list--\n");
    mexPrintf("       paico_mex(0)                 Converts the prefix tree to two arrays of patterns, where the index \n");
    mexPrintf("                                     stored with the prefix pattern is the index of the pattern in the array.\n");
    mexPrintf("                                     NOTE: This operation is irreversible.\n");
    mexPrintf(" [list1, list2] = paico_mex(index)  Returns the lists (patterns) associated to the index. \n");
    mexPrintf("                                     Assumes paico_mex(0) is run.\n");
    mexPrintf("\n");
    mexPrintf("Copyright 2011, Sami Hanhijärvi <sami.hanhijarvi@iki.fi>. See paico_mex.cpp for license.\n");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    static Node * tree = 0;
    static vector<vector<int> > highArray, lowArray;
    static int maxIndex = -1;   
    
    if (nrhs == 0 && nlhs == 0) {
        if (tree)
            delete tree;
        highArray.clear();
        lowArray.clear();
        
        tree = new Node(-1);
        return;
    } else if (nrhs == 0 && nlhs == 1) {
        string s;
        if (tree)
            s = tree->toString();
        plhs[0] = mxCreateString(s.c_str());
    } else if (nrhs == 1) {
        int index = (int)mxGetScalar(prhs[0]);
        if (index <= 0) {
            // Construct an array list of the prefix tree
            lowArray.resize(maxIndex+1);
            highArray.resize(maxIndex+1);
            convert(tree, lowArray, highArray);
        } else if (nlhs == 2) {            
            unsigned int len = lowArray[index].size();
            plhs[0] = mxCreateNumericMatrix(1,len,mxINT32_CLASS,mxREAL);            
            int *data = (int *)mxGetPr(plhs[0]);
            for (unsigned int i = 0; i < len; ++i) {
                // Matlab starts indexing from 1
                data[i] = lowArray[index][i]+1;
            }
            
            len = highArray[index].size();
            plhs[1] = mxCreateNumericMatrix(1,len,mxINT32_CLASS,mxREAL);            
            data = (int *)mxGetPr(plhs[1]);
            for (unsigned int i = 0; i < len; ++i) {
                // Matlab starts indexing from 1
                data[i] = highArray[index][i]+1;
            }
        } else {
            printUsage();
            mexErrMsgTxt("Illegal arguments. Check usage from above.");             
        }
            
    } else if (nrhs == 2) {
        double *input1 = mxGetPr(prhs[0]);
        int len1 = mxGetNumberOfElements(prhs[0]);
        
        double *input2 = mxGetPr(prhs[1]);
        int len2 = mxGetNumberOfElements(prhs[1]);
        
        if (nlhs == 0) {
            // Store the list
            if (!tree)
                tree = new Node(-1);
            unsigned int index = (unsigned int)input2[0];
            tree->store(input1, len1, index);
            if ((maxIndex < 0) || index >= maxIndex)
                maxIndex = index;
        } else if (nlhs == 1) {
            // Compare the lists
            int index = -1;
            if (tree)
                index = tree->findCommon(input1, len1, input2, len2);
            plhs[0] = mxCreateDoubleScalar(index);
        } else 
            mexErrMsgTxt("Illegal arguments. Check usage from above."); 
        
    } else {
        printUsage();
        mexErrMsgTxt("Illegal arguments. Check usage from above."); 
    }     
}
