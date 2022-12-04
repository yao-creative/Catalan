"""Author: Yi Yao Tan"""
import numpy as np
from typing import List, Tuple, Dict, Union, Any
import time
import matplotlib.pyplot as plt
#Notes I don't know if the copy is actually useful
#-----------------------------------------------

class Catalan_Structure:
    """Generic randomly generated combinatorial object counted by Catalan numbers of class size n"""
    def __init__(self,n: int):
        """n is the number of nodes in the planar tree"""
        self.n: int = n
        self.dyck_path: DyckPath = DyckPath([])
        self.plane_tree: Tree = Tree(0, [])
        self.binary_tree: Binary_Tree = Binary_Tree(0)
        self.non_crossing_partition: Noncrossing= Noncrossing([])
        self.heights: List[List[Tuple[int,int]]] = []
        self.plane_tree_graph: TreeGraph = TreeGraph()
        self.binary_tree_graph: TreeGraph = TreeGraph()
        self.times: Dict[str, float] = {"dyck":0,"noncrossing":0,"plane_tree":0, "plane_tree_graph":0, "binary_tree_graph":0}
        

    def generate(self, method: int= 1):
        """Actually generate the dyck path, binary tree, and plane tree"""
        start = time.time()
        self.dyck_path: DyckPath = self.generate_dyck_path(method) #has to be called first.
        stop = time.time()
        self.times["dyck"] = stop - start
        
        start = time.time()
        self.heights, self.non_crossing_partition = self.to_heights_and_noncrossing()
        stop = time.time()
        self.times["noncrossing"] = stop - start
        
        start  = time.time()
        self.plane_tree: Tree = self.to_planetree()
        stop = time.time()
        self.times["plane_tree"] = stop - start
        
        start  = time.time()
        self.plane_tree_graph: TreeGraph = self.to_planetree_graph()
        stop = time.time()
        self.times["plane_tree_graph"] = stop - start
        
        # start  = time.time()
        # self.binary_tree: Binary_Tree = self.to_binary_tree()
        # stop = time.time()
        # self.times["binary_tree"] = stop - start
        
        start  = time.time()
        self.binary_tree_graph: TreeGraph = self.to_binary_tree_graph()
        stop = time.time()
        self.times["binary_tree_graph"] = stop - start
        
        
        
    def show_times(self):
        print(f"Times\n {self.times}")
        
    def to_mathematica(self):
        self.dyck_path_mathematica: str = self.dyck_path.to_mathematica()
        self.non_crossing_partition_mathematica: str = self.non_crossing_partition.to_mathematica()
        self.plane_tree_mathematica: str = self.plane_tree.to_mathematica()
        # self.binary_tree_mathematica: str = self.binary_tree.to_mathematica()
        self.plane_tree_graph_mathematica: str = self.plane_tree_graph.to_mathematica()
        self.binary_tree_graph_mathematica: str = self.binary_tree_graph.to_mathematica()
        
    def show(self):
        print(f"Structures Generated:\nself.dyck_path {self.dyck_path}\nself.heights: {self.heights}\nself.non_crosing_partition: {self.non_crossing_partition}\nself.plane_tree: {self.plane_tree}\nself.plane_tree_graph: {self.plane_tree_graph} \nself.binary_tree_graph: {self.binary_tree_graph}")
        
    def show_mathematica(self):
        print(f"self.dyck_path_mathematica: {self.dyck_path_mathematica}\nself.non_crossing_partition_mathematica: {self.non_crossing_partition_mathematica}\nself.plane_tree_mathematica: {self.plane_tree_mathematica}\n \nself.plane_tree_graph_mathematica: {self.plane_tree_graph_mathematica} \nself.binary_tree_graph_mathematica: {self.binary_tree_graph_mathematica}")
        

        
    def random_permutation(self,n: int):
        return np.random.permutation(n)

    def random_permutation2(self,n: int):
        tab = [i for i in range(n)]
        for i in range(1,n):
            j = np.random.randint(0,i)
            tab[i], tab[j] = tab[j], tab[i]
        return np.array(tab)
    
    def generate_dyck_path(self, method: int = 1):
        """ generate a random path dyck path of length 2n + 1"""
        if method ==1:
            sigma = self.random_permutation(2*self.n +1)
        else:
            sigma = self.random_permutation2(2*self.n +1)
        w = []
        v = []
        u = []
        for i in range(2*self.n+1):
            if sigma[i] < self.n :
                w.append(1)
            else:
                w.append(0)
        curr_height, min, ind = 0, 0, 0
        for i in range(2*self.n + 1):
            if w[i] == 1:
                curr_height += 1
            else:
                curr_height -= 1
                if curr_height < min:
                    min = curr_height
                    ind = i
       
        u = w[ind+1:] + w[:ind+1]
        out = DyckPath(u)
        return out

    def check(self):
        """w is a path of length 2n+1 that should have n+1 down steps and n up steps""" 
        print(f"count 0: {self.dyck_path.path.count(0)} count 1: {self.dyck_path.path.count(1)}")
    
    def to_heights_and_noncrossing(self):
        """convert a path to a set of height classes of glued nodes using created path to represent the path, also generates non-crossing partition"""
        heights = []
        stack = []
        curr_height = 0
        noncrossing: List[Tuple[int,int]]  = []
        for i in range(len(self.dyck_path.path)-1): #ignore last one, just 2n in the beginning
            if self.dyck_path.path[i] == 1:
                curr_height += 1
                if curr_height > len(heights) :
                    heights.append([])
                stack.append(i)
            else:
                curr_height -= 1
                heights[curr_height].append((stack.pop(),i))
                noncrossing.append(heights[curr_height][-1])
                

        noncrossing_out = Noncrossing(noncrossing)
        return heights, noncrossing_out
    
    def to_planetree_graph(self):
        """converts a path to a plane tree in TreeGraph format"""
        stack  = [0]
        counter = 1
        TG = TreeGraph()
        for i in self.dyck_path.path:
            # print(f"i: {i} stack: {stack} TreeGraph: {TG}")
            if i == 1:
                TG.add_edge(stack[-1],counter)
                stack.append(counter)
                counter += 1
            if i == 0:
                stack.pop()
        return TG
            

    def to_binary_tree(self):
        """convert a heights representation to a binary tree of size n, dfs left to right label order"""
        heights = self.heights
        root = Binary_Tree((-1,2*self.n))
        if heights == []:
            return root #empty tree
        root.set_left(self.to_binary_tree_rec(heights,heights[0][0])) #first node planted without right child
        return root
    
    def to_binary_tree_rec(self, heights: list, tup: tuple):
        """ function: recursive helper to create a binary tree from a path
            heights: list of lists of tuples where each outer index denotes the height 
                    and each inner list stores the tuples of glued nodes
            tup: parent tuple of glued node
            label: label of the node"""
        root = Binary_Tree(tup)
        
        if heights == []: #top right corner
            root.set_left(Binary_Leaf("left"))
            root.set_right(Binary_Leaf("right"))
            return root

        else:
            #cases for left child
            if len(heights) > 1:
                new_heights = heights[1:]
                for i in range(len(new_heights[0])): 
                    #find the left child which is an inner vertex of the tree.
                    # 1) find in the next height level a child node belonging to the parent within the tuple range
                    # 2) has the lowest index
                    if new_heights[0][i][0] > tup[0] and new_heights[0][i][1] < tup[1]:
                        root.set_left(self.to_binary_tree_rec(new_heights, heights[1][i]))
                        break
                #some cases without left child, a leaf or outer vertext
                if root.left == None:
                    root.set_left(Binary_Leaf("left"))
            else:
                root.set_left(Binary_Leaf("left")) #no greater height
                
            #cases for right child, if there are nodes to the right of the same height 
            if len(heights[0]) > 1:
                new_heights = [heights[0][1:]] + heights[1:]
                root.set_right(self.to_binary_tree_rec(new_heights, heights[0][1]))
            else: #since increasing every node of the same height are child of right of the parent.
                root.set_right(Binary_Leaf(1))
                

            return root
    def to_binary_tree_graph(self):
        """converts a dyck path to a binary tree graph in TreeGraph format"""    
        BTG = TreeGraph()
        stack = [0]
        counter = 1
        for i in self.dyck_path.path[:-1]:  
            # print(f"i: {i} stack: {stack}")
            if i == 1:
                BTG.add_edge(stack[-1],counter) #left edge
                stack.append(counter)
                counter += 1
            if i == 0:
                stack.pop()
                BTG.add_edge(stack[-1],counter) #add right edge
                stack[-1] = counter #replace top of stack with new highest right edge
                counter +=1
        return BTG
                
    def to_planetree(self):
        """convert a path to a planar tree by using a height class representation"""
        heights = self.heights
        root = Tree((-1,2*self.n)) #plant plane tree with root of the whole path
        children = []
        for i in range(len(heights[0])):
            children.append(self.to_planetree_rec((heights[1:]), (heights[0][i])))  
        root.set_children(children)
        return root
    
    def to_planetree_rec(self, heights: list, tup): 
        """ function: recursive helper to create a planar tree from a path
            heights: list of lists of tuples where each outer index denotes the height 
                    and each inner list stores the tuples of glued nodes
            tup: parent tuple of glued node"""
          
        if heights == []:
            root = Leaf(tup)    
            return root
        else:    
            children = []
            for i in range(len(heights[0])):
                if heights[0][i][0] > tup[0] and heights[0][i][1] < tup[1]: 
                    #this is a child. Parent is (a,b) and child is (c,d) if a < c < d < b, 
                    #thus generate children subtrees with next new height level, new parent is the child
                    children.append(self.to_planetree_rec((heights[1:]), (heights[0][i])))
                    
            if len(children)> 0:
                #set the children for the Tree
                root = Tree(tup)
                root.set_children(children)
                return root
            else:
                root = Leaf(tup)
                return root


"""Trees can have children both trees and leaves, and leaves have no children"""
#------------------Generic Tree Class------------------
class Node:
    def __init__(self, value, children = []):
        """value is a tuple of glued nodes (start, end)
            children is a list of subtrees"""
        self.value = value
        self.children = children

class TreeGraph:
    def __init__(self):
        self.edge_list = []
        self.num_nodes = 0
        
    def __str__(self):
        return f"TreeGraph({self.edge_list})"
    
    def add_edge(self, edge1:int, edge2:int):
        self.edge_list.append((edge1, edge2))
        self.num_nodes +=1
        
    def to_mathematica(self):
        edge_list_mathematica = []
        for edge in self.edge_list:
            edge_list_mathematica.append(f"{edge[0]} -> {edge[1]}")
        nodes_list = ",".join([str(i) for i in range(self.num_nodes)])
        return f"TreeGraph ["+ "{" + ",".join(edge_list_mathematica) + '}, GraphLayout -> \"LayeredDigraphEmbedding\"]'
        
class Tree(Node):
    def __init__(self, value, children = []):
        super().__init__(value, children)
        
    def __str__(self):
        children_str = ",".join((str(child) for child in self.children))
        return f"Tree(%s, [%s])" % (self.value, children_str)
    
    def to_mathematica(self):
        children_str = ",".join([child.to_mathematica() for child in self.children])
        # print(f"children_str: {children_str}")
        return f"Tree[a, {{%s}}]" % (children_str)
        
    def set_children(self, new_children: list):
        """set list of children trees"""
        self.children: list = new_children
        
    def add_child(self, child = None):
        self.children.append(child)    

class Leaf(Node):
    def __init__(self, value):
        super().__init__(value)
        
    def __str__(self):
        return f"Leaf(%s)" % (str(self.value))
    
    def to_mathematica(self):
        return f"a" # % (str(self.value))

##------------------Binary Tree Class------------------
class Binary_Node:
    def __init__(self, value, left = None, right = None):
        """value is a tuple of glued nodes (start, end)
            left is a subtree
            right is a subtree"""
        self.value = value
        self.left = left
        self.right = right
        
    def __str__(self):
        return f"Binary_Node(%s, %s, %s)" % (str(self.value), str(self.left), str(self.right))
    
    def to_mathematica(self):
        return "Binary_Node[%s, %s, %s]" % (str(self.value), str(self.left), str(self.right))
    
class Binary_Tree(Binary_Node):
    def __init__(self, value, left = None, right = None):
        super().__init__(value, left, right)
        
    def __str__(self):
        return f"Binary_Tree(%s, %s, %s)" % (str(self.value), str(self.left), str(self.right))
    
    def to_mathematica(self):
        if self.left is not None:
            left = self.left.to_mathematica()
            children = ""
            if self.right is not None: 
                right = self.right.to_mathematica() 
                children = f"{{{left}, {right}}}"
            else: 
                children = f"{{{left}}}"
            return f"Tree[ , {children}]"
        return f"Tree[ , {{}}]"
        
    def set_left(self, left: Binary_Node):
        self.left = left
        
    def set_right(self, right: Binary_Node):
        self.right = right
    
class Binary_Leaf(Binary_Node):
    def __init__(self, value):
        super().__init__(value)
        
    def __str__(self):
        return f"Binary_Leaf(%s)" % (str(self.value))
    
    def to_mathematica(self):
        return f" " #% (str(self.value))
class Noncrossing:
    def __init__(self, edge_list: List[Tuple[int,int]] ):
        self.edge_list = edge_list
        self.num_nodes = len(edge_list)*2
        
    def __str__(self):
        return f"Noncrossing({self.edge_list})"
    
    def to_mathematica(self):
        edge_str = ", ".join([f"UndirectedEdge[{edge[0]}, {edge[1]}]" for edge in self.edge_list])
        vertex_str = ", ".join([str(i) for i in range(self.num_nodes)])
        return f"Graph[{{{vertex_str}}}, {{{edge_str}}}, GraphLayout -> {{\"CircularEmbedding\", \"OptimalOrder\" -> False}}, VertexLabels -> \"Name\"]"

class DyckPath:
    def __init__(self, path: List[int], ys: List[int] =[]) :
        self.path = path
        self.ys = self.generate_ys() if len(ys) == 0 else ys
        
    def generate_ys(self) -> List[int]:
        counter =0
        ys = [0]
        for i in self.path:
            if i == 1:
                counter +=1
                ys.append(counter)
            else:
                counter -=1
                ys.append(counter)
        return ys
                
    def __str__(self):
        return f"DyckPath({self.path})"
    
    def to_mathematica(self):
        return f"{{{','.join([str(i) for i in self.path])}}}"
    
    def plot(self):
        plt.title("Dyck Path")
        print(f"self.ys: {self.ys}")
        plt.xlim(0,len(self.ys)-2)
        plt.ylim(0,max(self.ys)+1)
        plt.plot(self.ys)
        plt.show()
    
if __name__ == "__main__":
    n = 5
    A = Catalan_Structure(n)
    A.generate(2)
    # A.dyck_path.plot()
    print(f"A.dyck_path.ys: {A.dyck_path.ys} length: {len(A.dyck_path.ys)}")
    A.to_mathematica()
    A.show_mathematica()

    
        