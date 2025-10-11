# MoonBit Graph Theory Library (DS)

This library provides a comprehensive collection of graph data structures and algorithms implemented in MoonBit, supporting various graph analysis and operations.

## Features

- **Graph Structures**: Support for directed and undirected graphs, with weighted or unweighted edges
- **Basic Operations**: Add/remove nodes and edges, access neighbor nodes
- **Graph Traversal**: Depth-First Search (DFS), Breadth-First Search (BFS), Topological Sort
- **Path Algorithms**: Dijkstra, Bellman-Ford, Floyd-Warshall, A*
- **Connectivity Analysis**: Connected Components, Strongly Connected Components, Bridges and Cut Vertices
- **Tree Algorithms**: Minimum Spanning Tree (Kruskal, Prim), Tree Diameter
- **Flow Networks**: Maximum Flow, Minimum Cost Maximum Flow, Minimum Feasible Flow
- **Feature Detection**: Cycle Detection, Bipartite Graph Detection, Euler Path/Circuit
- **Graph Metrics**: Degree Centrality, PageRank, Clustering Coefficient, Graph Density
- **Graph Coloring**: Greedy Coloring, DSatur Coloring

## Getting Started

### Creating a Graph and Adding Nodes and Edges

```moonbit
    // Create an undirected graph
    let g : Graph2[String, Int] = Graph2::new(directed=false)

    // Add nodes
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")

    // Add weighted edges
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=2) |> ignore
    g.add_edge(a, c, weight=3) |> ignore

    // Get basic graph information
    println(g.node_count())  // 3
    println(g.edge_count())  // 3
```

### Graph Traversal

#### DFS (Depth-First Search)

```moonbit
    // Perform DFS traversal
    let dfs = Dfs::new(g, start_node)
    let mut node = dfs.next(g)
    while node != None {
    let current = node.unwrap()
    // Process current node
    node = dfs.next(g)
    }
```

#### BFS (Breadth-First Search)

```moonbit
    // Perform BFS traversal
    let bfs = Bfs::new(g, start_node)
    let mut node = bfs.next(g)
    while node != None {
    let current = node.unwrap()
    // Process current node
    node = bfs.next(g)
    }
```

#### Topological Sort

```moonbit
    // Perform topological sorting
    let topo = Topo::new(g)
    let mut node = topo.next(g)
    while node != None {
    let current = node.unwrap()
    // Process nodes in topological order
    node = topo.next(g)
    }
```

### Shortest Path Algorithms

#### Dijkstra's Algorithm

```moonbit
    // Find single-source shortest paths
    let distances = g.dijkstra(start_node)
    for (node, distance) in distances {
    println("Distance to node \{node.ix}: \{distance}")
    }
```

#### Bellman-Ford Algorithm (Supports Negative Weights)

```moonbit
    // Find single-source shortest paths, supporting negative edge weights
    let (distances, has_negative_cycle) = g.bellman_ford(start_node)
    if has_negative_cycle {
    println("Graph contains a negative cycle")
    } else {
    // Process distances array
    }
```

#### A* Pathfinding Algorithm

```moonbit
    // Define a heuristic function
    fn heuristic(a : NodeIndex, b : NodeIndex) -> Int {
    // Calculate heuristic distance estimate
    // ...
    return estimated_distance
    }

    // Find path from start to goal
    let path = g.a_star(start, goal, heuristic)

    // Process the path array
    for node in path {
    // Visit each node in the path
    }
```

### Connectivity Analysis

#### Connected Components (Undirected Graph)

```moonbit
    let components = g.connected_components()
    println("Number of connected components: \{components.length()}")

    // View nodes in each component
    for (i, component) in components.enumerate() {
    println("Component \{i} contains \{component.length()} nodes")
    }
```

#### Strongly Connected Components (Directed Graph)

```moonbit
    let sccs = g.strongly_connected_components()
    println("Number of strongly connected components: \{sccs.length()}")

    // Process each SCC
    for scc in sccs {
    // Process nodes in each strongly connected component
    }
```

#### Finding Bridges and Cut Vertices

```moonbit
    // Find all bridges
    let bridges = g.find_bridges()

    // Find all cut vertices (articulation points)
    let articulation_points = g.find_articulation_points()
```

### Minimum Spanning Tree

#### Kruskal's Algorithm

```moonbit
    // Get the MST edge set
    let mst = g.kruskal_mst()

    // Calculate total MST weight
    let total_weight = mst.fold(0, fn(acc, edge) { acc + edge.2 })
```

#### Prim's Algorithm

```moonbit
    // Compute MST using Prim's algorithm
    let mst = g.prim_mst()

    // Process MST edges
    for (u, v, weight) in mst {
    // Process each edge in the MST
    }
```

### Flow Network Algorithms

#### Maximum Flow

```moonbit
    // Calculate maximum flow from source to sink
    let max_flow = g.max_flow(source, sink)
```

#### Minimum Cost Maximum Flow

```moonbit
    // Define edge costs
    let costs = @hashmap.new()
    costs.set((node1, node2), cost)

    // Calculate minimum cost maximum flow
    let (flow, cost) = g.min_cost_max_flow(source, sink, costs)
```

### Graph Feature Analysis

#### Cycle Detection

```moonbit
    // Find all cycles in the graph
    let cycles = g.find_cycles()

    // Check each cycle
    for cycle in cycles {
    // Process each cycle
    }
```

#### Bipartite Graph Detection

```moonbit
    // Check if graph is bipartite
    let (is_bipartite, colors) = g.is_bipartite()
    if is_bipartite {
    // Get the two parts of the bipartite graph
    let left = Array::new()
    let right = Array::new()
    for i in 0..<g.node_count() {
        if colors[i] == 0 {
        left.push(i)
        } else if colors[i] == 1 {
        right.push(i)
        }
    }
    }
```

#### Euler Path/Circuit

```moonbit
    // Check for Euler circuit (all vertices have even degree)
    let has_circuit = g.has_euler_circuit()

    // Check for Euler path (at most two vertices with odd degree)
    let has_path = g.has_euler_path()
```

### Graph Metrics

#### Degree Centrality

```moonbit
    // Calculate degree centrality for each node
    let centrality = g.degree_centrality()

    // Find the node with highest centrality
    let mut max_centrality = 0.0
    let mut max_node = -1
    for i in 0..<centrality.length() {
    if centrality[i] > max_centrality {
        max_centrality = centrality[i]
        max_node = i
    }
    }
```

#### PageRank

```moonbit
    // Calculate PageRank
    let ranks = g.pagerank(0.85, 100, 0.0001)

    // Process PageRank values for each node
    for i in 0..<ranks.length() {
    println("PageRank of node \{i}: \{ranks[i]}")
    }
```

#### Clustering Coefficient

```moonbit
    // Calculate the average clustering coefficient of the graph
    let cc = g.clustering_coefficient()
```

### Graph Coloring

```moonbit
    // Color the graph using greedy algorithm
    let (color_count, colors) = g.greedy_coloring()

    // Color the graph using DSatur algorithm
    let (color_count, colors) = g.dsatur_coloring()
```

## Advanced Usage

### Compressed Graph Representation

```moonbit
    // Create compressed graph representation to optimize memory usage
    let compressed = g.to_compressed()

    // Use the compressed graph
    let neighbors = compressed.get_neighbors(node)
    ```

    ### Graph Diameter

    ```moonbit
    // Calculate tree diameter (longest path)
    let (diameter, path) = g.tree_diameter()
```

## Notes

- Nodes and edges can have weights of various types
- After deleting nodes, node indices may change and need careful handling
- Use `|> ignore` to handle function calls that don't need return values

## Advanced Examples

### Complex Network Analysis

```moonbit
    // Create a complex network
    let g = Graph2::new(directed=true)

    // Add nodes and edges
    // ...

    // Compute strongly connected components
    let sccs = g.strongly_connected_components()

    // Calculate PageRank for each SCC
    for scc in sccs {
    let sub_graph = g.subgraph(scc)
    let ranks = sub_graph.pagerank(0.85, 100, 0.0001)
    // Analyze PageRank results for each subgraph
    }
```
