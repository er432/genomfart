from functools import partial
from collections import deque
import networkx as nx
import numpy as np
import itertools
import heapq

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
            else:
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
    def run_LBP_rbp(self, maxit=1000):
        """ Runs loopy belief propagation using residual belief propagation
        
        maxit: int
            Maximum number of iterations to perform
        """
        # Reset messages
        self.reset_nonobs_messages()
        ## Send initial messages up and down using BFS
        ## Also, create queue of (-message_priority, message_edge, message_name)
        # Initialize all priority to 0
        message_queue = []
        ## Create dictionary of (message_edge, message_name) -> current priority
        priority_dict = {}
        root_node = sorted(self.nodes(),key=self.degree)[-1]
        bfs = list(nx.bfs_edges(self,source=root_node))
        for to_node, from_node in reversed(bfs):
            if self.node[from_node]['bipartite'] == 'v':
                self.send_message_to_factor(from_node,to_node)
                heapq.heappush(message_queue,(0.,(from_node,to_node),'to_f_message'))
                priority_dict[((from_node,to_node),'to_f_message')] = -0.001
            else:
                self.send_message_to_var(from_node, to_node)
                heapq.heappush(message_queue,(0.,(from_node,to_node),'to_v_message'))
                priority_dict[((from_node,to_node),'to_v_message')] = -0.001
        for from_node, to_node in bfs:
            if self.node[from_node]['bipartite'] == 'v':
                self.send_message_to_factor(from_node,to_node)
                heapq.heappush(message_queue,(0.,(from_node,to_node),'to_f_message'))
                priority_dict[((from_node,to_node),'to_f_message')] = -0.001
            else:
                self.send_message_to_var(from_node, to_node)
                heapq.heappush(message_queue,(0.,(from_node,to_node),'to_v_message'))
                priority_dict[((from_node,to_node),'to_v_message')] = -0.001
        # Iterate based on priority
        for it in xrange(maxit):
            # Get highest priority edge
            neg_priority,(n1,n2),message_name = heapq.heappop(message_queue)
            current_message_val = np.array(self.edge[n1][n2][message_name])
            # Send the new message
            var_node = n1 if self.node[n1]['bipartite'] == 'v' else n2
            factor_node = n1 if self.node[n1]['bipartite'] == 'f' else n2
            to_node = None
            from_node = None
            if message_name == 'to_f_message':
                self.send_message_to_factor(var_node,factor_node)
                to_node = factor_node
                from_node = var_node
            else:
                self.send_message_to_var(factor_node,var_node)
                to_node = var_node
                from_node = factor_node
            # Calculate new node priority
            priority = np.max(np.abs(np.log(np.divide(self.edge[n1][n2][message_name],
                                                   current_message_val))))
            if np.isnan(priority):
                priority=0.
            # Put into queue
            heapq.heappush(message_queue,(-priority,(n1,n2),message_name))
            priority_dict[((from_node,to_node),message_name)] = -priority
            # Update neighboring nodes, not including the from node
            neighbor_message_name = 'to_v_message' if message_name == 'to_f_message' else 'to_f_message'
            for n in self.neighbors(to_node):
                if n == to_node: continue
                # Remove the message from the queue
                message_queue.remove((priority_dict[((to_node,n),neighbor_message_name)],
                                          (to_node,n),neighbor_message_name))
                # Update message
                current_message_val = np.array(self.edge[to_node][n][neighbor_message_name])
                if neighbor_message_name == 'to_f_message':
                    self.send_message_to_factor(to_node, n)
                else:
                    self.send_message_to_var(to_node,n)
                # Calculate residual
                priority = np.max(np.abs(np.log(np.divide(self.edge[to_node][n][neighbor_message_name],
                                                              current_message_val))))
                if np.isnan(priority):
                    priority=0.
                # Put back into queue
                heapq.heappush(message_queue,(-priority,(to_node,n),neighbor_message_name))
                priority_dict[((to_node,n),neighbor_message_name)]=-priority
            # End if all priorities are 0
            if all([x[0]==0. for x in message_queue]):
                break
    def send_message_to_factor(self, var, factor_name):
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
            self.edge[var][factor_name]['f_message_sent'] = True
            return
        message = np.ones_like(self.edge[var][factor_name]['to_f_message'])
        # Multiply out the message from factors to compute the message
        for n in self.neighbors(var):
            if n != factor_name:
                message = np.multiply(message,
                                      self.edge[n][var]['to_v_message'])
        # Send the message
        self.edge[var][factor_name]['to_f_message'] = message
        self.edge[var][factor_name]['f_message_sent'] = True
    def send_message_to_var(self, factor_name, var):
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
            self.edge[factor_name][var]['v_message_sent'] = True
            return
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
        # Send the message
        self.edge[factor_name][var]['to_v_message'] = message
        self.edge[factor_name][var]['v_message_sent'] = True
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
    def convert_to_discrete_factor_graph(self):
        """ Converts the existing bayesian network to a factor graph
        
        Returns
        -------
        A factor_graph object
        """
        factor = discrete_factor_graph()
        def factor_func(node_cpd, *args):
            return node_cpd.get_prob(node_value=args[-1], 
                    ordered_arg_vals=tuple(args[:-1]))
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
        return factor

