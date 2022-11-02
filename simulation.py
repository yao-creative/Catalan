"""Author: Yi Yao Tan"""
from codecs import namereplace_errors
import numpy as np
import copy

#Notes I don't know if the copy is actually useful
#-----------------------------------------------
class Catalan_Structure:
    """Generic randomly generated combinatorial object counted by Catalan numbers of class size n"""
    def __init__(self,n: int):
        """n is the number of nodes in the planar tree"""
        self.n: int = n
        self.dyck_path: list = self.generate_dyck_path() #has to be called first.
        self.plane_tree: Tree = self.to_planetree()
        self.binary_tree: Binary_Tree = self.to_binary_tree()
        print(f"self.dyck_path {self.dyck_path} \nself.plane_tree: {self.plane_tree} \nself.binary_tree: {self.binary_tree}")

    def random_permutation(self,n: int):
        return np.random.permutation(n)
    
    def generate_dyck_path(self):
        """ generate a random path dyck path of length 2n + 1"""
        sigma = self.random_permutation(2*self.n +1)
        w = []
        u = []
        for i in range(2*self.n+1):
            if sigma[i] < self.n :
                w.append(1)
            else:
                w.append(0)
        # print(f"generated w count 0: {w.count(0)} count 1: {w.count(1)}")
        curr_height, min, ind = 0, 0, 0
        for i in range(2*self.n + 1):
            if w[i] == 1:
                curr_height += 1
            else:
                curr_height = curr_height - 1
                if curr_height < min:
                    min = curr_height
                    ind = i
                
        print(f"w: {w} min: {min} ind: {ind}")
        #everything after ind is positive walk since it is minimum
        u = w[ind+1:] + w[:ind+1]
        return u

    def check(self, w: list):
        """w is a path of length 2n+1 that should have n+1 down steps and n up steps""" 
        print(f"count 0: {w.count(0)} count 1: {w.count(1)}")
    
    def to_heights(self):
        """convert a path to a set of height classes of glued nodes using created path to represent the path"""
        heights = []
        stack = []
        curr_height = 0
        for i in range(len(self.dyck_path)-1): #ignore last one, just 2n in the beginning
            if self.dyck_path[i] == 1:
                curr_height += 1
                if curr_height > len(heights) :
                    heights.append([])
                stack.append(i)
            else:
                curr_height -= 1
                heights[curr_height].append((stack.pop(),i))

        return heights
    
    def to_binary_tree(self):
        """convert a heights representation to a binary tree of size n, dfs left to right label order"""
        heights = self.to_heights()
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
            root.set_left(Binary_Leaf(0))
            root.set_right(Binary_Leaf(1))
            return root

        else:
            #cases for left child
            if len(heights) > 1:
                new_heights = heights[1:]
                for i in range(len(new_heights[0])): #find the left child:
                    # 1) belonging to the parent within the tuple range
                    # 2) has the lowest index
                    if new_heights[0][i][0] > tup[0] and new_heights[0][i][1] < tup[1]:
                        root.set_left(self.to_binary_tree_rec(new_heights, heights[1][i]))
                        break
                #some cases without left child:
                if root.left == None:
                    root.set_left(Binary_Leaf(0))
            else:
                root.set_left(Binary_Leaf(0)) #no greater height
                
            #cases for right child, if there are nodes to the right of the same height 
            if len(heights[0]) > 1:
                new_heights = [heights[0][1:]] + heights[1:]
                root.set_right(self.to_binary_tree_rec(new_heights, heights[0][1]))
            else: #since increasing every node of the same height are child of right of the parent.
                root.set_right(Binary_Leaf(1))
                

            return root
                
        
    def to_planetree(self):
        """convert a path to a planar tree by using a height class representation"""
        heights = self.to_heights()
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
                if heights[0][i][0] > tup[0] and heights[0][i][1] < tup[1]: #this is a child
                    children.append(self.to_planetree_rec((heights[1:]), (heights[0][i])))
            if len(children)> 0:
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
    
    
        
class Tree(Node):
    def __init__(self, value, children = []):
        super().__init__(value, children)
        
    def __str__(self):
        children_str = ",".join((str(child) for child in self.children))
        return f"Tree(%s, [%s])" % (self.value, children_str)
    
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
    
class Binary_Tree(Binary_Node):
    def __init__(self, value, left = None, right = None):
        super().__init__(value, left, right)
        
    def __str__(self):
        return f"Binary_Tree(%s, %s, %s)" % (str(self.value), str(self.left), str(self.right))
    
    def set_left(self, left: Binary_Node):
        self.left = left
        
    def set_right(self, right: Binary_Node):
        self.right = right
    
class Binary_Leaf(Binary_Node):
    def __init__(self, value):
        super().__init__(value)
        
    def __str__(self):
        return f"Binary_Leaf(%s)" % (str(self.value))
    



if __name__ == "__main__":
    n = 50
    Catalan_Structure(n)
    # print(f"tree child 0: {tree.children[0]}")

    
        