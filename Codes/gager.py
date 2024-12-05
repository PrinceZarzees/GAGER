import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import networkx as nx
import sys

# Replace these file paths with the actual paths to your data
# take input from command line
expression_before = sys.argv[1]
expression_after = sys.argv[2]


# Read gene expression matrices
before = pd.read_csv(expression_before, sep=',',index_col=0)
after = pd.read_csv(expression_after, sep=',',index_col=0)
before=before.T
after=after.T

# Function to calculate linear regression
def linear_regression(x, a, b):
    return a * x + b


# Function to perform regression analysis
def perform_regression(gene1, gene2,regression_function):
    # Take average expression for each gene across matrices
    if gene1 not in before.index:
      avg_gene1_before=0
    else:
      avg_gene1_before = before.loc[gene1].mean()
    if gene1 not in after.index:
      avg_gene1_after=0
    else:
      avg_gene1_after = after.loc[gene1].mean()

    if gene2 not in before.index:
      avg_gene2_before=0
    else:
      avg_gene2_before = before.loc[gene2].mean()
    if gene2 not in after.index:
      avg_gene2_after=0
    else:
      avg_gene2_after = after.loc[gene2].mean()


    # Calculate the relative expression changes
    x = [avg_gene1_before,avg_gene1_after]
    y = [avg_gene2_before,avg_gene2_after]
    params, covariance = curve_fit(regression_function, x, y)

    return params



# Read data from CSV files
file1_path = sys.argv[3]  # Replace with the actual file path for the first graph
file2_path = sys.argv[4]  # Replace with the actual file path for the second graph
expression1_path = sys.argv[1]  # Replace with the actual file path for the first gene expression matrix
expression2_path = sys.argv[2]  # Replace with the actual file path for the second gene expression matrix

# Read CSV files into Pandas DataFrames
graph1_df = pd.read_csv(file1_path)
graph2_df = pd.read_csv(file2_path)

# Read gene expression matrices into Pandas DataFrames

file1 = pd.read_csv(expression1_path, sep=',', index_col=0)
file2 = pd.read_csv(expression2_path, sep=',', index_col=0)

file1=file1.T
file2=file2.T


all_genes = list(set(file1.index).union(file2.index))


file1 = file1.reindex(all_genes, fill_value=0)
file2 = file2.reindex(all_genes, fill_value=0)

# Initialize an empty dataframe for the result
result = pd.DataFrame(index=all_genes)

# Iterate through each gene
for gene in result.index:
    # Calculate average expression values for each gene in both files
    avg_expr_file1 = file1.loc[gene].mean()
    avg_expr_file2 = file2.loc[gene].mean()

    # Find the expression difference
    expr_difference = avg_expr_file2 - avg_expr_file1

    # Calculate standard deviation for each gene in both files
    std_expr_file1 = file1.loc[gene].std()
    std_expr_file2 = file2.loc[gene].std()

    # Add the expression difference to the result dataframe
    result.at[gene, 'Expression Difference'] = expr_difference



# Sort edges by importance and take the top 10%
# Define the threshold for edge importance
# Normalize the edge weights for graph1
graph1_df['normalized_importance'] = graph1_df['importance'] / max(graph1_df['importance'])
graph2_df['normalized_importance'] = graph2_df['importance'] / max(graph2_df['importance'])
# print (graph1_df)

threshold = 0.10

# Correct the column names for source and target nodes
important_edges_graph1 = graph1_df[graph1_df['normalized_importance'] > threshold][['TF', 'target', 'normalized_importance']]
important_edges_graph2 = graph2_df[graph2_df['normalized_importance'] > threshold][['TF', 'target', 'normalized_importance']]

# Create NetworkX graphs from the filtered edges with weights
graph1 = nx.from_pandas_edgelist(important_edges_graph1, source='TF', target='target', edge_attr='normalized_importance')
graph2 = nx.from_pandas_edgelist(important_edges_graph2, source='TF', target='target', edge_attr='normalized_importance')

# Get the common edges
common_edges = set(graph1.edges).intersection(graph2.edges)

# Remove common edges from both graphs
graph1.remove_edges_from(common_edges)
graph2.remove_edges_from(common_edges)

# Create the difference graph by composing the filtered graphs
# difference_graph = nx.compose(graph1, graph2)
difference_graph=graph2

all_nodes = set(graph1_df['TF']).union(set(graph1_df['target'])).union(set(graph2_df['TF'])).union(set(graph2_df['target']))
# difference_graph.add_nodes_from(all_nodes)

# Compute the topological sorting
diff_graph = nx.DiGraph(difference_graph.edges())
diff_graph.add_nodes_from(difference_graph.nodes())
# Define a custom sorting key
def custom_sort_key(node):
    # Get the outdegree of the node in the difference graph
    total_weight=0
    for child in diff_graph.successors(node):
      total_weight+=difference_graph[node][child]['normalized_importance']
    # Return a tuple with outdegree as the first element and total_weight as the second element
    return -total_weight

def remove_cycles(graph):
    """
    Removes cycles from a directed graph by discarding edges with the lowest importance in each cycle.

    Parameters:
    graph (nx.DiGraph): A directed graph with an 'importance' attribute for edges.

    Returns:
    None: Modifies the graph in place to remove cycles.
    """
    while not nx.is_directed_acyclic_graph(graph):
        # Find all simple cycles in the graph
        cycles = list(nx.simple_cycles(graph))
        
        if not cycles:
            break  # Exit if no cycles are found (shouldn't happen since the graph is not a DAG)
        
        # Process each cycle
        for cycle in cycles:
            # Get all edges in the current cycle
            edges_in_cycle = [(cycle[i], cycle[(i + 1) % len(cycle)]) for i in range(len(cycle))]
            
            # Find the edge with the lowest importance in the cycle
            min_edge = min(
                edges_in_cycle,
                key=lambda edge: graph[edge[0]][edge[1]].get('normalized_importance', float('inf'))
            )
            
            # Remove the edge with the lowest importance
            graph.remove_edge(*min_edge)
            break  # Recompute cycles after modifying the graph to avoid stale cycles

    # Check if the graph is acyclic after processing
    if nx.is_directed_acyclic_graph(graph):
        print("Cycles removed successfully. The graph is now a DAG.")
    else:
        print("Failed to remove cycles completely.")


remove_cycles(diff_graph)
print (diff_graph)




# topologically_sorted_nodes = list(nx.topological_sort(diff_graph.subgraph(largest_component)))
topologically_sorted_nodes = list(nx.topological_sort(diff_graph))


# Initialize a dictionary to store the layer information for each node
layer_info = {node: 0 for node in topologically_sorted_nodes}

# Iterate over the topologically sorted nodes
for node in topologically_sorted_nodes:
    # Get the predecessors (parents) of the current node
    parents = list(diff_graph.predecessors(node))
    # Assign the layer based on the number of parents
    if len(parents) == 0:
        # If the node has no parents, it belongs to layer 0
        layer_info[node] = 0
    else:
        # If the node has one or more parents, it belongs to the next layer
        max_parent_layer = max(layer_info[parent] for parent in parents)
        layer_info[node] = max_parent_layer + 1

# Define a custom sorting function to sort nodes within each layer
def sort_nodes_within_layer(layer_nodes):
    return sorted(layer_nodes, key=custom_sort_key)

# Sort the nodes within each layer
for layer in set(layer_info.values()):
    # Get the nodes belonging to the current layer
    layer_nodes = [node for node, layer_value in layer_info.items() if layer_value == layer]

    # Sort the nodes within the current layer
    sorted_layer_nodes = sort_nodes_within_layer(layer_nodes)

    # Update the topologically sorted nodes list with the sorted nodes within the current layer
    for i, node in enumerate(sorted_layer_nodes):
        topologically_sorted_nodes.remove(node)
        topologically_sorted_nodes.insert(i + layer * len(sorted_layer_nodes), node)

def find_downstream_nodes(node, difference_graph):
    """Find the downstream nodes affected by a change in the given node."""
    downstream_nodes = set()
    for successor in diff_graph.successors(node):
        downstream_nodes.add(successor)
        downstream_nodes.update(find_downstream_nodes(successor, diff_graph))
    return downstream_nodes
def calculate_total_difference(result):
    """Calculate the total node expression difference for all downstream nodes."""
    total_difference = 0
    for node in topologically_sorted_nodes:
      try:
        # if abs(result.loc[node,'Expression Difference'])>result.loc[node,'sigma']:
          total_difference += abs(result.loc[node,'Expression Difference'])
      except:
        pass
    return total_difference
g_result=result.copy()
# Initialize variables to keep track of the best node and its total expression difference
result=g_result.copy()
TFs_to_change = []
best_total_difference = calculate_total_difference(result)

# Initialize threshold for stopping criteria
print (best_total_difference)


threshold = best_total_difference*0.35   
  # Adjust as needed
def parents_of(node, graph):
    parents = []
    for edge in graph.in_edges(node):
        parents.append(edge[0])
    return parents

# Iterate through the topologically sorted nodes
for node in topologically_sorted_nodes:
    # Make the expression of the current node 0
    # For example, if node is a key in the expression dictionary, set its value to 0
    temp = result.copy()
    try:
      if abs(temp.loc[node,'Expression Difference'])<0.5*temp.loc[node,'sigma']:
        continue
    except:
      pass

    # Find the downstream nodes affected by this change
    downstream_nodes = find_downstream_nodes(node, difference_graph)
 
    # Calculate the expression difference for each downstream node
    temp.at[node, 'Expression Difference'] = 1e-5
    for d_node in downstream_nodes:

      parent_list=parents_of(d_node,diff_graph)

      sum_ax=0
      for parent_node in parent_list:
        params = perform_regression(parent_node, d_node,linear_regression)
    
        sum_ax+=params[0]*temp.at[parent_node,'Expression Difference']*difference_graph[parent_node][d_node]['normalized_importance']
        
      temp.at[d_node, 'Expression Difference'] = (sum_ax)/len(parent_list)

      # params = perform_regression(node, d_node,linear_regression)
      # temp.at[d_node, 'Expression Difference'] = 0



    # Calculate the total node expression difference for all downstream nodes
    total_difference = calculate_total_difference(temp)

    # Check if the total node expression difference is less than the previous value
    if total_difference < best_total_difference:
        # Update the previous value and save the node
        best_total_difference = total_difference
        TFs_to_change.append(node)
        result=temp.copy()

    # Check if the stopping criterion is met
    if total_difference < threshold:
        break

# Print the best node and its total expression difference
print(TFs_to_change)
print(len(TFs_to_change))
print("Total expression difference:", best_total_difference)

