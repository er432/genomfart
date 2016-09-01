from multiprocessing.dummy import Pool
from functools import partial
from collections import deque
from numba import jit
import networkx as nx
import numpy as np
import itertools
import heapq

def write_libpgm_json(json_dict,write_file):
    num_func = lambda x: ("%0.3g" % x)[1:] if 1.e-12< x<1 else ("%0.3g" % np.abs(np.round(x)))
    write_file.write('{\n')
    write_file.write("\t\"V\": [%s],\n" % ", ".join(map(lambda x: "\"%s\"" % x, 
                        json_dict["V"])))
    for i,e in enumerate(json_dict['E']):
        write_line = None
        n_edges = len(json_dict['E'])
        if i == 0:
            if n_edges > 1:
                write_line = "\t\"E\": [[%s],\n" % ", ".join(map(lambda x: "\"%s\"" % x, e))
            else:
                write_line = "\t\"E\": [[%s]],\n" % ", ".join(map(lambda x: "\"%s\"" % x, e))
        else:
            if i == (n_edges-1):
                write_line = "\t\t[%s]],\n" % ", ".join(map(lambda x: "\"%s\"" % x, e))
            else:
                write_line = "\t\t[%s],\n" % ", ".join(map(lambda x: "\"%s\"" % x, e))
        write_file.write(write_line)
    write_file.write('\t"Vdata": {\n')
    n_nodes = len(json_dict['Vdata'])
    for i,(k,v) in enumerate(json_dict['Vdata'].items()):
        write_file.write('\t\t"%s": {\n' % k)
        write_file.write('\t\t\t"numoutcomes": %d,\n' % v["numoutcomes"])
        write_file.write('\t\t\t"vals": [%s],\n' % \
                             ", ".join(map(lambda x: "\"%s\"" % x, v['vals'])))
        if v['parents'] is None:
            write_file.write('\t\t\t"parents": None,\n')
        else:
            write_file.write('\t\t\t"parents": [%s],\n' % \
                             ", ".join(map(lambda x: "\"%s\"" % x, v['parents'])))
        if v['children'] is None:
            write_file.write('\t\t\t"children": None,\n')
        else:
            write_file.write('\t\t\t"children": [%s],\n' % \
                             ", ".join(map(lambda x: "\"%s\"" % x, v['children'])))
        if type(v['cprob']) == list:
            write_file.write('\t\t\t"cprob": [%s]\n' % \
                                 ", ".join(map(num_func, v['cprob'])))
        else:
            n_probs = len(v['cprob'])
            for j,(k1,v1) in enumerate(v['cprob'].items()):
                if j == 0:
                    write_file.write('\t\t\t"cprob": {\n')
                    write_line = "\t\t\t\t\"%s\": [%s],\n" % \
                          (k1,
                               ", ".join(map(num_func, v1)))
                else:
                    if j == (n_probs-1):
                        write_line = "\t\t\t\t\"%s\": [%s]\n" % \
                              (k1,
                                   ", ".join(map(num_func, v1)))
                    else:
                        write_line = "\t\t\t\t\"%s\": [%s],\n" % \
                              (k1,
                                   ", ".join(map(num_func, v1)))
                write_file.write(write_line)
            write_file.write("\t\t\t}\n")
        if i == (n_nodes-1):
            write_file.write("\t\t}\n")
        else:
            write_file.write("\t\t},\n")
    write_file.write("\t}\n")
    write_file.write("}\n")

## Function for sending messages in parallel
def _par_send_message(G,from_node,to_node):
    if G.node[to_node]['bipartite'] == 'v' and \
      G.node[to_node]['observed']: return None
    if G.node[from_node]['bipartite'] == 'v':
        return (((from_node,to_node),'to_f_message'),
                    G.send_message_to_factor(from_node,to_node,update=False))
    else:
        return (((from_node,to_node),'to_v_message'),
                    (G.send_message_to_var(from_node, to_node,update=False)))

class Messenger(object):
    def __init__(self, G):
        self._G = G
    def __call__(self,e):
        return _par_send_message(self._G,e[0],e[1])

class cpd(object):
    """ Base class for a conditional probability distribution
    """
    def __init__(self, node_name, parents):
        """ Instantiates the CPD
        
        Parameters
        ----------
        node_name : hashable
            Name of the node
        parents : tuple or list
            Names of the parents
        """
        self._node_name = node_name
        self._parents = tuple(parents)
        ## Create dictionary of parent -> index
        self._parent_ind_dict = {k:i for i,k in enumerate(parents)}
    def get_name(self):
        """ Gets the name of the node
        
        Returns
        -------
        Name of the node
        """
        return self._node_name
    def get_parents(self):
        """ Gets the name of the node's parents
        
        Returns
        -------
        Names of the parents
        """
        return self._parents

class discrete_cpd(cpd):
    """ Class used to represent a conditional probability
    distribution with only discrete values
    """
    def __init__(self, node_name, parents, parent_value_dict, node_values,
                     prob_dict):
        """ Instantiates the discrete CPD
        
        Parameters
        ----------
        node_name : hashable
            Name of the node
        parents : tuple or list
            Ordered names of the parents
        parent_value_dict : dict
            Dictionary of parent -> (Ordered values parent can take)
        node_values : tuple
            Ordered values this child node can take
        prob_dict : dict
            Dictionary of (parent1 value, parent2 value, ...) ->  
            np.array([prob child value1, prob child value2, ...])
        """
        super(discrete_cpd,self).__init__(node_name, parents)
        ## Error checking ##b
        # Check that number of parents conforms to prob_dict keys
        if not all(map(lambda k: len(k) == len(self._parents), prob_dict.keys())):
            raise ValueError("Number of parents not conforming to prob_dict key length")
        # Check that number of node values conforms to prob_dict vector length
        if not all(map(lambda v: len(v) == len(node_values),prob_dict.values())):
            raise ValueError("Length of probability vector not conforming to node_values length")
        # Enforce sum-to-1 constraint
        prob_dict = {k:np.array(v)/np.sum(v) for k,v in prob_dict.items()}
        self._parent_value_dict = parent_value_dict
        self._node_values = node_values
        ## Create dictionary of node_value -> index
        self._node_ind_dict = {k:i for i,k in enumerate(node_values)}
        self._prob_dict = prob_dict
    def get_prob(self,ordered_arg_vals=None, log=False, node_value=None, **parent_vals):
        """ Gets the (log) probability of node values given the values
        of its parents
        
        Parameters
        ----------
        ordered_arg_vals : tuple
            Ordered argument values for parents
        log : bool
            Whether to return the log-probability rather than the actual probability
        node_value : None or value
            If not None, a specified node value at which the CPD should be evaluated
        parent_vals :
            keyword arguments of parent_name -> value
            
        Returns
        -------
        Conditional (log) probability. If node_value is None, will return the whole
        probability vector. Otherwise, gives the probability for the given node_value
        entry
        """
        if ordered_arg_vals is not None:
            parent_vals = {self._parents[i]:v for i,v in enumerate(ordered_arg_vals)}
        if len(parent_vals) != len(self._parents):
            raise ValueError("The number of parent values given does not equal the number of \
            parents for this CPD")
        # Get the key
        key = np.zeros(len(self._parents)).tolist()
        for k,v in parent_vals.items():
            key[self._parent_ind_dict[k]] = v
        key = tuple(key)
        if node_value is not None:
            if log:
                return np.log(self._prob_dict[key][self._node_ind_dict[node_value]])
            else:
                return self._prob_dict[key][self._node_ind_dict[node_value]]
        else:
            if log:
                return np.log(self._prob_dict[key])
            else:
                return self._prob_dict[key]

class sparse_discrete_cpd(discrete_cpd):
    """ Use to represent and discrete CPD where most elements are zero
    """
    def __init__(self, node_name, parents, parent_value_dict, node_values,
                     prob_dict):
        """ Instantiates the discrete CPD
        
        Parameters
        ----------
        node_name : hashable
            Name of the node
        parents : tuple or list
            Ordered names of the parents
        parent_value_dict : dict
            Dictionary of parent -> (Ordered values parent can take)
        node_values : tuple
            Ordered values this child node can take
        prob_dict : dict
            Dictionary of (parent1 value, parent2 value, ...) ->  
            np.array([prob child value1, prob child value2, ...])
        """
        super(sparse_discrete_cpd,self).__init__(node_name, parents,
                      parent_value_dict, node_values, prob_dict)
        ## Create dictionary of parent_vals -> plausible prob inds
        self._plausibility_dict = {}
        for parent_vals, child_probs in self._prob_dict.iteritems():
            plausible_inds = []
            for i,prob in enumerate(child_probs):
                if prob > 1.e-12:
                    plausible_inds.append(i)
            self._plausibility_dict[parent_vals] = tuple(plausible_inds)
    def get_plausible_val_iter(self):
        """ Get a generator of dictionary of node_name -> value for each
        present combination

        Returns
        -------
        Generator of dictionary of {parent or child name -> value}
        """
        for parent_vals, child_probs in self._prob_dict.iteritems():
            base_dict = dict(zip(self._parents, parent_vals))
            for i in self._plausibility_dict[parent_vals]:
                base_dict[self._node_name] = self._node_values[i]
                yield dict(base_dict)
    def get_prob(self,ordered_arg_vals=None, log=False, node_value=None, **parent_vals):
        """ Gets the (log) probability of node values given the values
        of its parents
        
        Parameters
        ----------
        ordered_arg_vals : tuple
            Ordered argument values for parents
        log : bool
            Whether to return the log-probability rather than the actual probability
        node_value : None or value
            If not None, a specified node value at which the CPD should be evaluated
        parent_vals :
            keyword arguments of parent_name -> value
            
        Returns
        -------
        Conditional (log) probability. If node_value is None, will return the whole
        probability vector. Otherwise, gives the probability for the given node_value
        entry
        """
        if ordered_arg_vals is not None:
            parent_vals = {self._parents[i]:v for i,v in enumerate(ordered_arg_vals)}
        if len(parent_vals) != len(self._parents):
            raise ValueError("The number of parent values given does not equal the number of \
            parents for this CPD")
        # Get the key
        key = np.zeros(len(self._parents)).tolist()
        for k,v in parent_vals.items():
            key[self._parent_ind_dict[k]] = v
        key = tuple(key)
        if key not in self._prob_dict:
            if log:
                return -np.inf
            else:
                return 0.
        if node_value is not None:
            if log:
                return np.log(self._prob_dict[key][self._node_ind_dict[node_value]])
            else:
                return self._prob_dict[key][self._node_ind_dict[node_value]]
        else:
            if log:
                return np.log(self._prob_dict[key])
            else:
                return self._prob_dict[key]

class discrete_factor_graph(nx.Graph):
    """ Class used to represent a factor graph for discrete variables, which is a bipartite graph
    with factor and variable nodes. Node attributes under bipartite are 
    'f' and 'v', respectively. Factor nodes have attribute 'func', containing
    their functions
    """
    def __init__(self):
        """ Instantiates the factor graph
        """
        super(discrete_factor_graph, self).__init__()
    def add_factor_node(self, factor_name, factor_vars, factor_func,
                            var_value_dict):
        """ Adds a factor node to the graph, along with its variable nodes

        Parameters
        ----------
        factor_name : hashable
            Name to give the factor node
        factor_vars : tuple or list
            Ordered list of variables connected to the factor node. This should
            be the same order as the arguments given to the factor_func
        factor_func : function
            Function taking ordered factor_vars values as arguments and returning
            some scalar value
        var_value_dict : dict
            Dictionary of variable -> (Ordered values variable can take)
        """
        ## Create factor node
        super(discrete_factor_graph, self).add_node(factor_name, bipartite='f',
           func=factor_func, 
           arg_order={k:i for i,k in enumerate(factor_vars)},
           ordered_args = factor_vars,
           var_value_dict=var_value_dict)
        ## Connect to variable nodes
        for var in factor_vars:
            # Add the variable node
            super(discrete_factor_graph,self).add_node(var, bipartite='v', 
                  state=None,
                  observed=False,
                  val_order_dict = {k:i for i,k in enumerate(var_value_dict[var])})
            # Add edge with message passing attributes
            super(discrete_factor_graph,self).add_edge(factor_name,var,
              to_v_message=np.ones(len(var_value_dict[var])),
              to_f_message=np.ones(len(var_value_dict[var])),
              v_message_sent = False,
              f_message_sent = False)
    def get_norm_beliefs(self, var):
        """ Gets the normalized beliefs for a variable node based on the current
        messages running to it
        
        Parameters
        ----------
        var : hashable
            The variable to get the beliefs for
        
        Returns
        -------
        Vector of normalized beliefs (probabilities)
        """
        # Get product of messages
        message_prod = np.ones(len(self.node[var]['val_order_dict']))
        for n in self.neighbors(var):
            message_prod = np.multiply(message_prod, 
                           self.edge[n][var]['to_v_message'])
        return (message_prod/np.sum(message_prod))
    def get_random_spanning_tree(self):
        """ Gets a random spanning tree of this graph

        Returns
        -------
        A random spanning tree of this factor graph
        """
        tree = nx.Graph()
        # Choose a random starting node
        current_node = np.random.choice(self.nodes())
        while len(tree) < len(self):
            # Randomly select a neighbor
            neighbor = np.random.choice(self.neighbors(current_node))
            # Add the edge from the current node to the neighbor if the
            # neighbor hasn't been encountered before
            if neighbor not in tree:
                tree.add_edge(current_node, neighbor)
            current_node = neighbor
        return tree
    def reset_nonobs_messages(self):
        """ Resets all messages on the graph that aren't involving observed
        variables
        """
        # Ensure all messages at uniform
        for n1,n2 in self.edges_iter():
            if self.node[n1]['bipartite'] == 'v':
                if self.node[n1]['observed']: continue
            elif self.node[n2]['bipartite'] == 'v':
                if self.node[n2]['observed']: continue
            self.edge[n1][n2]['to_v_message'] = np.ones_like(self.edge[n1][n2]['to_v_message']) 
            self.edge[n1][n2]['to_f_message'] = np.ones_like(self.edge[n1][n2]['to_f_message']) 
            self.edge[n1][n2]['v_message_sent'] = False
            self.edge[n1][n2]['f_message_sent'] = False
    def run_BP_tree(self):
        """ Run belief propagation for a tree. Will throw an error if not
        a tree

        Raises
        ------
        ValueError 
            If not a tree
        """
        # Ensure graph is a tree
        if not nx.is_tree(self):
            raise ValueError("Factor graph is not a tree")
        # Reset messages
        self.reset_nonobs_messages()
        # Create message queue
        message_queue = deque(sorted(self.nodes(), 
                                                 key=lambda n: self.degree(n)))
        # Send messages
        while len(message_queue) > 0:
            # Get a node
            node = message_queue.popleft()
            # Get the node type
            node_type = self.node[node]['bipartite']
            to_node_message_type = ("%s_message_sent" % node_type)
            from_node_message_type = ("%s_message_sent" % ('v' if node_type == 'f' else 'f'))
            # Get neighbors that haven't sent message to this node
            unsent_from_message_nodes = frozenset([n for n in self.neighbors(node) if \
                          self.edge[n][node][to_node_message_type] == False])
            # Get neighbors that haven't been sent messages from this node
            unsent_to_message_nodes = [n for n in self.neighbors(node) if \
                            self.edge[node][n][from_node_message_type] == False]
            n_to_send = len(unsent_to_message_nodes)
            # Get node degree
            node_degree = self.degree(node)
            # Tell this node to send message to unsent nodes if it has been sent
            # messages from all other nodes
            if len(unsent_from_message_nodes) <= 1:
                for send_node in unsent_to_message_nodes:
                    if all((len(unsent_from_message_nodes) == 1, 
                                send_node in unsent_from_message_nodes)):
                        # The send node is the one that hasn't sent the message
                        pass
                    elif len(unsent_from_message_nodes) == 1:
                        # The send node is NOT the one that hasn't sent the message
                        continue
                    if node_type == 'v':
                        self.send_message_to_factor(node, send_node)
                    else:
                        self.send_message_to_var(node, send_node)
                    n_to_send -= 1
            # Check if this node still has messages to send
            if n_to_send > 0:
                message_queue.append(node)
    def run_synchronous_LBP(self, maxit=100, verbose=False):
        """ Runs loopy belief propagation by synchronously updating
        messages while they change

        Parameters
        ----------
        maxit : int
            Maximum number of iterations to perform
        """ 
        # Run an initial synchronous update on all edges
        self.run_synchronous_update()
        # Run additional updates
        to_update_list = []
        for i in xrange(maxit):
            changed = None
            if i == 0:
                changed = self.run_synchronous_update(return_changed=True)
            else:
                changed = self.run_synchronous_update(edge_yielder=to_update_list,
                                                          return_changed=True)
            # If nothing has changed, break from loop
            if len(changed) == 0:
                break
            to_update_list = []
            added_edges = set()
            # Add all messages emanating from nodes that have received changed
            # messages to the update list
            for from_node,to_node in changed:
                for n in self.neighbors(to_node):
                    if (to_node,n) not in added_edges:
                        to_update_list.append((to_node,n))
                        added_edges.add((to_node,n))
        if verbose:
            if i == (maxit-1):
                print("Did not completely converge after %d iterations" % (maxit))
            else:
                print("Converged after %d iterations" % i)
    def run_synchronous_update(self, edge_yielder=None, return_changed = False):
        """ Updates all nodes simultaneously
        
        Parameters
        ----------
        edge_yielder : iterable
            Lists edges to update
        return_changed : bool
            Whether to return list of edges that were modified
        """ 
        # Create dictionary of ((from_node,to_node),message_type) -> message
        message_dict = {}
        run_all = (edge_yielder is None)
        if edge_yielder is None:
            edge_yielder = self.edges_iter()
        message_dict.update(**dict(filter(lambda x: x is not None,
                                map(Messenger(self),
                           edge_yielder))))
        if run_all:
            message_dict.update(**dict(filter(lambda x: x is not None,
                              map(Messenger(self),
                           map(lambda x: tuple(reversed(x)),self.edges_iter())))))
        # Send all messages, and potentially add to the list of messages
        # that were changed
        return_list = []
        for ((n1,n2),message_type),message in message_dict.items():
            if return_changed:
                if not np.allclose(self.edge[n1][n2][message_type],message):
                    return_list.append((n1,n2))
            self.edge[n1][n2][message_type] = message
        if return_changed:
            return return_list
    def send_message_to_factor(self, var, factor_name, update=True):
        """ Sends a message from a variable to a factor
        
        Parameters
        ----------
        var : hashable
            Variable name
        factor_name : hashable
            Factor name
        """
        if not all((self.node[var]['bipartite'] == 'v', 
                        self.node[factor_name]['bipartite'] == 'f')):
            raise ValueError("Nodes are not a variable and a factor")
        # Skip if observed (message should already be set)
        if self.node[var]['observed']:
            if update:
                self.edge[var][factor_name]['f_message_sent'] = True
                return
            else:
                return self.edge[factor_name][var]['to_f_message']
        message = np.ones_like(self.edge[var][factor_name]['to_f_message'])
        # Multiply out the message from factors to compute the message
        for n in self.neighbors(var):
            if n != factor_name:
                message = np.multiply(message,
                                      self.edge[n][var]['to_v_message'])
        # Ensure moderate values
        message /= np.sum(message)
        if 0 < np.min(message) < 1.e-10:
            message /= np.max(message)
        elif np.max(message) > 1.e10:
            message /= np.max(message)
        if update:
            # Send the message
            self.edge[var][factor_name]['to_f_message'] = message
            self.edge[var][factor_name]['f_message_sent'] = True
        else:
            return message
    def send_message_to_var(self, factor_name, var, update=True):
        """ Sends a message from a factor to a variable
        
        Parameters
        ----------
        var : hashable
            Variable name
        factor_name : hashable
            Factor name
        """
        if not all((self.node[var]['bipartite'] == 'v', 
                        self.node[factor_name]['bipartite'] == 'f')):
            raise ValueError("Nodes are not a variable and a factor")
        # Skip if variable is observed
        if self.node[var]['observed']:
            if update:
                self.edge[factor_name][var]['v_message_sent'] = True
                return
            else:
                return self.edge[factor_name][var]['to_v_message']
        # Instantiate the message block
        dims = np.zeros(self.degree(factor_name),dtype=np.int64)
        dim_iter = []
        for n in self.neighbors(factor_name):
            dims[self.node[factor_name]['arg_order'][n]] = \
              len(self.edge[n][factor_name]['to_f_message'])
            dim_iter.append(tuple(np.arange(len(self.edge[n][factor_name]['to_f_message'])).tolist()))
        dim_iter = itertools.product(*dim_iter)
        message = np.ones(dims)
        # Fill in the values, moving through all combos of values
        for var_args in \
          itertools.product(*[self.node[factor_name]['var_value_dict'][k] for k in \
                                  self.node[factor_name]['ordered_args']]):
            # Get the indices for this var combo
            inds = next(dim_iter)
            # Multiply by the factor function value
            message[inds] *= self.node[factor_name]['func'](*var_args)
            # Multiply through with non-recipient arg values
            for i,n in enumerate(self.node[factor_name]['ordered_args']):
                if n != var:
                    message[inds]*=self.edge[n][factor_name]['to_f_message'][inds[i]]
        # Marginalize over non-recipient variables
        sub_num = 0
        for i,n in enumerate(self.node[factor_name]['ordered_args']):
            if n != var:
                message = np.sum(message, axis=(i-sub_num))
                sub_num += 1
        # Ensure moderate values
        message /= np.sum(message)
        if 0 < np.min(message) < 1.e-10:
            message /= np.max(message)
        elif np.max(message) > 1.e10:
            message /= np.max(message)
        if update:
            # Send the message
            self.edge[factor_name][var]['to_v_message'] = message
            self.edge[factor_name][var]['v_message_sent'] = True
        else:
            return message
    def set_variable_state(self, var, val):
        """ Sets a variable's state, given observed data
        
        Parameters
        ----------
        var : hashable
            The variable name
        val : hashable
            The observed value
       """
        if not self.node[var]['bipartite'] == 'v':
            raise ValueError("%s is not a variable node" % str(var))
        self.node[var]['state'] = val
        self.node[var]['observed'] = True
        # Set all messages to reflect observed status
        for n in self.neighbors(var):
            self.edge[var][n]['to_f_message'] = np.zeros_like(self.edge[var][n]['to_f_message'])
            self.edge[var][n]['to_f_message'][self.node[var]['val_order_dict'][val]] = 1
            self.edge[var][n]['f_message_sent'] = True
    def set_variable_unobserved(self, var):
        """ Sets a variable's state to unobserved

        Parameters
        ----------
        var : hashable
            The variable name
        """
        if not self.node[var]['bipartite'] == 'v':
            raise ValueError("%s is not a variable node" % str(var))
        self.node[var]['state'] = None
        self.node[var]['observed'] = False
        for n in self.neighbors(var):
            self.edge[var][n]['to_f_message'] = np.ones_like(self.edge[var][n]['to_f_message'])
            self.edge[var][n]['f_message_sent'] = False

class sparse_discrete_factor_graph(discrete_factor_graph):
    """ Class used to represent a factor graph for discrete variables if probabilities are sparse, 
    It is a bipartite graph
    with factor and variable nodes. Node attributes under bipartite are 
    'f' and 'v', respectively. Factor nodes have attribute 'func', containing
    their functions
    """
    def __init__(self):
        """ Instantiates the factor graph
        """
        super(sparse_discrete_factor_graph, self).__init__()
    def send_message_to_var(self, factor_name, var, update=True):
        """ Sends a message from a factor to a variable
        
        Parameters
        ----------
        var : hashable
            Variable name
        factor_name : hashable
            Factor name
        """
        if not all((self.node[var]['bipartite'] == 'v', 
                        self.node[factor_name]['bipartite'] == 'f')):
            raise ValueError("Nodes are not a variable and a factor")
        # Skip if variable is observed
        if self.node[var]['observed']:
            if update:
                self.edge[factor_name][var]['v_message_sent'] = True
                return
            else:
                return self.edge[factor_name][var]['to_v_message']
        # Instantiate the message as zeros
        message = np.zeros_like(self.edge[factor_name][var]['to_v_message'])
        # Fill in the values, moving through all combos of values that are non-zero
        for var_arg_dict in self.node[factor_name]['val_iter_func']():
            # Get ordered arguments
            var_args = map(lambda k: var_arg_dict[k], self.node[factor_name]['ordered_args'])
            # Start with the conditional probability
            message_val_part = self.node[factor_name]['func'](*var_args)
            if message_val_part < 1.e-12: 
                 #Skip if the probability is zero to start with
                continue
            # Go through, multiplying by the appropriate parts of the 
            # to factor messages from non-recipient args
            for i,n in enumerate(self.node[factor_name]['ordered_args']):
                if n != var:
                    var_ind = self.node[n]['val_order_dict'][var_arg_dict[n]]
                    mult_factor = self.edge[n][factor_name]['to_f_message'][var_ind]
                    if mult_factor > 1.e-12:
                        message_val_part *= mult_factor
                    else:
                        message_val_part = 0.
                        break
            # Add message part to message
            message[self.node[var]['val_order_dict'][var_arg_dict[var]]] += message_val_part
        # Ensure moderate values
        message /= np.sum(message)
        if 0 < np.min(message) < 1.e-10:
            message /= np.max(message)
        elif np.max(message) > 1.e10:
            message /= np.max(message)
        if update:
            # Send the message
            self.edge[factor_name][var]['to_v_message'] = message
            self.edge[factor_name][var]['v_message_sent'] = True
        else:
            return message

class bayes_net(nx.DiGraph):
    """ Class used to represent a Bayesian network
    """
    def __init__(self):
        """ Instantiates the Bayesian network
        """
        super(bayes_net, self).__init__()
    def add_node(self, cpd):
        """ Adds a node to the Bayes net
        
        Parameters
        ----------
        cpd : cpd object
            A conditional probability distribution for the node
        """
        # Add the base node, placing the CPD object there
        super(bayes_net, self).add_node(cpd.get_name(), cpd=cpd,
                                            state=None, observed=False)
        # Add the edges
        for parent in cpd.get_parents():
            super(bayes_net, self).add_edge(parent, cpd.get_name())
    def convert_to_discrete_factor_graph(self, sparse=False):
        """ Converts the existing bayesian network to a factor graph
        
        Parameters
        ----------
        sparse : bool
            Whether the CPDs in the factor graph should be regarded as sparse

        Returns
        -------
        A factor_graph object
        """
        factor = None
        if sparse:
            factor = sparse_discrete_factor_graph()
        else:
            factor = discrete_factor_graph()
        def factor_func(node_cpd, *args):
            return node_cpd.get_prob(node_value=args[-1], 
                    ordered_arg_vals=tuple(args[:-1]))
        def val_iter_func(node_cpd, *args):
            return node_cpd.get_plausible_val_iter()
        # Go through nodes, creating factor for each node
        for node in self.nodes_iter():
            # Get ordered node parents from CPD
            ordered_parents = self.node[node]['cpd'].get_parents()
            # Create a node name of (parents, node)
            node_name = tuple(list(ordered_parents)+[node])
            factor_vars = node_name
            var_value_dict = dict(self.node[node]['cpd']._parent_value_dict)
            var_value_dict[node] = self.node[node]['cpd']._node_values
            # Create appropriate factor function
            factor.add_factor_node(node_name, factor_vars,
                   partial(factor_func,self.node[node]['cpd']),
                   dict(var_value_dict))
            # If sparse, add function that will generate plausible values
            if sparse:
                factor.node[node_name]['val_iter_func'] = partial(val_iter_func,
                                                        self.node[node]['cpd'])
        return factor
    def write_libpgm_file(self, filename):
        """ Writes to a file that can be read by libpgm

        filename : str
            File to write
        """
        json_dict = {"V":[],"E":[],"Vdata":{}}
        for node in self.nodes_iter():
            json_dict["V"].append(node)
            json_dict["Vdata"][node] = {}
            json_dict["Vdata"][node]["numoutcomes"] = len(self.node[node]['cpd']._node_values)
            json_dict["Vdata"][node]["vals"] = list(self.node[node]['cpd']._node_values)
            if len(self.predecessors(node)) == 0:
                json_dict["Vdata"][node]["parents"] = None
            else:
                json_dict["Vdata"][node]["parents"] = list(self.node[node]['cpd'].get_parents())
            if len(self.successors(node)) == 0:
                json_dict["Vdata"][node]["children"] = None
            else:
                json_dict["Vdata"][node]["children"] = self.successors(node)
            if len(self.node[node]['cpd']._prob_dict) == 1:
                json_dict["Vdata"][node]["cprob"] = self.node[node]['cpd']._prob_dict.values()[0].tolist()
            else:
                json_dict["Vdata"][node]["cprob"] = {}
                for k,v in self.node[node]['cpd']._prob_dict.items():
                    json_dict["Vdata"][node]["cprob"][str(list(k))] = v.tolist()
        for e1,e2 in self.edges_iter():
            json_dict["E"].append([e1,e2])
        handle = open(filename,'w')
        write_libpgm_json(json_dict,handle)
        handle.close()
        return
                
            

