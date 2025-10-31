# MoonBit Graph Theory Library (DS)

This library provides a comprehensive collection of graph data structures and algorithms implemented in MoonBit, supporting various graph analysis and operations.

## Features

- **Graph Structures**: Support for directed and undirected graphs, with weighted or unweighted edges
- **Basic Operations**: Add/remove nodes and edges, access neighbor nodes, graph density analysis
- **Graph Traversal**: Depth-First Search (DFS), Breadth-First Search (BFS), Topological Sort
- **Path Algorithms**: 
  - Single-source: Dijkstra, Bellman-Ford, A*
  - All-pairs: Floyd-Warshall, Johnson
  - Multiple paths: K shortest paths
- **Connectivity Analysis**: Connected Components, Strongly Connected Components, Bridges and Cut Vertices
- **Tree Algorithms**: Minimum Spanning Tree (Kruskal, Prim), Second MST, Tree Diameter
- **Flow Networks**: 
  - Maximum Flow: Ford-Fulkerson, Dinic, Push-Relabel
  - Cost-based: Minimum Cost Maximum Flow, Minimum Cost Flow
  - Constrained: Minimum Feasible Flow
- **Feature Detection**: 
  - Cycles: Cycle Detection, Euler Path/Circuit Finding
  - Graph Types: Bipartite Graph Detection, Graph Isomorphism
  - Cliques: Maximum Clique Detection
- **Graph Metrics**: 
  - Centrality: Degree Centrality, PageRank
  - Structure: Clustering Coefficient, Graph Density
  - Path Analysis: Minimum Bottleneck Path
- **Graph Coloring**: Greedy Coloring, DSatur Coloring
- **Advanced Operations**: Node removal, compressed graph representation

## Getting Started

### Creating a Graph and Adding Nodes and Edges

```moonbit
test {
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
}
```

### Graph Traversal

#### DFS (Depth-First Search)

```moonbit
test {
    // Create a graph
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=2) |> ignore
    g.add_edge(c, d, weight=3) |> ignore
    g.add_edge(a, d, weight=4) |> ignore
    
    // Perform DFS traversal starting from node A
    let start_node = a
    let dfs = Dfs::new(g, start_node)
    let mut node = dfs.next(g)
    while node != None {
        let _current = node.unwrap()
        // Process current node
        node = dfs.next(g)
    }
}
```

#### BFS (Breadth-First Search)

```moonbit
test {
    // Create a graph
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=2) |> ignore
    g.add_edge(c, d, weight=3) |> ignore
    g.add_edge(a, d, weight=4) |> ignore
    
    // Perform BFS traversal starting from node A
    let start_node = a
    let bfs = Bfs::new(g, start_node)
    let mut node = bfs.next(g)
    while node != None {
        let _current = node.unwrap()
        // Process current node
        node = bfs.next(g)
    }
}
```

#### Topological Sort

```moonbit
test {
    // Create a directed acyclic graph
    let g : Graph2[String, Int] = Graph2::new(directed=true)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(a, c, weight=2) |> ignore
    g.add_edge(b, d, weight=3) |> ignore
    g.add_edge(c, d, weight=4) |> ignore
    
    // Perform topological sorting
    let topo = Topo::new(g)
    let mut node = topo.next(g)
    while node != None {
        let _current = node.unwrap()
        // Process nodes in topological order
        node = topo.next(g)
    }
}
```

### Shortest Path Algorithms

#### Dijkstra's Algorithm

```moonbit
test {
    // Create a weighted graph
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=4) |> ignore
    g.add_edge(a, c, weight=2) |> ignore
    g.add_edge(b, c, weight=1) |> ignore
    g.add_edge(c, d, weight=3) |> ignore
    g.add_edge(b, d, weight=5) |> ignore
    
    // Find single-source shortest paths from node A
    let start_node = a
    let _distances = g.dijkstra(start_node)
    // Note: distances is a HashMap, iterate through it
    // for (node, distance) in distances {
    //     println("Distance to node \{node.ix}: \{distance}")
    // }
}
```

#### Bellman-Ford Algorithm (Supports Negative Weights)

```moonbit
test {
    // Create a graph with potential negative weights
    let g : Graph2[String, Int] = Graph2::new(directed=true)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=4) |> ignore
    g.add_edge(a, c, weight=2) |> ignore
    g.add_edge(b, c, weight=-1) |> ignore  // Negative weight
    g.add_edge(c, d, weight=3) |> ignore
    g.add_edge(b, d, weight=5) |> ignore
    
    // Find single-source shortest paths, supporting negative edge weights
    let start_node = a
    let (_distances, has_negative_cycle) = g.bellman_ford(start_node)
    if has_negative_cycle {
        println("Graph contains a negative cycle")
    } else {
        // Process distances array
        // Note: distances is a HashMap, iterate through it
        // for (node, distance) in distances {
        //     println("Distance to node \{node.ix}: \{distance}")
        // }
    }
}
```

#### A* Pathfinding Algorithm

```moonbit
test {
    // Create a graph for pathfinding
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=4) |> ignore
    g.add_edge(a, c, weight=2) |> ignore
    g.add_edge(b, c, weight=1) |> ignore
    g.add_edge(c, d, weight=3) |> ignore
    g.add_edge(b, d, weight=5) |> ignore
    
    // Define a heuristic function (e.g., Manhattan distance)
    fn heuristic(_a : NodeIndex, _b : NodeIndex) -> Int {
        // Simple heuristic: return 0 for demonstration
        // In practice, you might use actual coordinates
        return 0
    }

    // Find path from start to goal
    let start = a
    let goal = d
    let path = g.a_star(start, goal, heuristic)

    // Process the path array
    for node in path {
        // Visit each node in the path
        println("Path node: \{g.node_weights[node.ix]}")
    }
}
```

#### Floyd-Warshall Algorithm (All-Pairs Shortest Paths)

```moonbit
test {
    // Create a graph for all-pairs shortest paths
    let g : Graph2[String, Int] = Graph2::new(directed=true)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=2) |> ignore
    g.add_edge(a, c, weight=5) |> ignore
    g.add_edge(c, d, weight=1) |> ignore
    
    // Calculate all-pairs shortest paths
    let distances = g.floyd_warshall()
    
    // Access distance between any two nodes
    let dist_a_to_c = distances[a.ix][c.ix]
    let dist_b_to_d = distances[b.ix][d.ix]
    
    if dist_a_to_c != None {
        println("Distance from A to C: \{dist_a_to_c.unwrap()}")
    }
    if dist_b_to_d != None {
        println("Distance from B to D: \{dist_b_to_d.unwrap()}")
    }
}
```

#### Johnson Algorithm (All-Pairs Shortest Paths with Negative Weights)

```moonbit
test {
    // Create a graph that may have negative weights
    let g : Graph2[String, Int] = Graph2::new(directed=true)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    let e = g.add_node(weight="E")
    let f = g.add_node(weight="F")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=2) |> ignore
    g.add_edge(a, d, weight=5) |> ignore
    g.add_edge(b, f, weight=3) |> ignore
    g.add_edge(f, e, weight=1) |> ignore
    g.add_edge(e, d, weight=1) |> ignore
    
    // Calculate all-pairs shortest paths using Johnson's algorithm
    let distances = g.johnson_apsp()
    
    // Check distances between various nodes
    let dist_a_to_c = distances[a.ix][c.ix]
    let dist_a_to_d = distances[a.ix][d.ix]
    let dist_a_to_e = distances[a.ix][e.ix]
    
    if dist_a_to_c != None {
        println("Distance from A to C: \{dist_a_to_c.unwrap()}")
    }
    if dist_a_to_d != None {
        println("Distance from A to D: \{dist_a_to_d.unwrap()}")
    }
    if dist_a_to_e != None {
        println("Distance from A to E: \{dist_a_to_e.unwrap()}")
    }
}
```

#### K Shortest Paths

```moonbit
test {
    // Create a graph with multiple paths
    let g : Graph2[String, Int] = Graph2::new(directed=true)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    let e = g.add_node(weight="E")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, e, weight=2) |> ignore
    g.add_edge(a, c, weight=3) |> ignore
    g.add_edge(b, d, weight=1) |> ignore
    g.add_edge(c, d, weight=2) |> ignore
    g.add_edge(d, e, weight=1) |> ignore
    
    // Find K shortest paths from A to E
    let start = a
    let goal = e
    let k = 3
    let paths = g.k_shortest_paths(start, goal, k)
    
    // Process each path
    for i in 0..<paths.length() {
        let (weight, path) = paths[i]
        println("Path \{i + 1}: weight = \{weight}, length = \{path.length()}")
    }
}
```

### Connectivity Analysis

#### Connected Components (Undirected Graph)

```moonbit
test {
    // Create a graph with multiple components
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    let e = g.add_node(weight="E")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=2) |> ignore
    g.add_edge(d, e, weight=3) |> ignore
    // Note: no edge between c and d, creating two components
    
    let components = g.connected_components()
    println("Number of connected components: \{components.length()}")

    // View nodes in each component
    for i in 0..<components.length() {
        let component = components[i]
        println("Component \{i} contains \{component.length()} nodes")
    }
}
```

#### Strongly Connected Components (Directed Graph)

```moonbit
test {
    // Create a directed graph with cycles
    let g : Graph2[String, Int] = Graph2::new(directed=true)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=2) |> ignore
    g.add_edge(c, a, weight=3) |> ignore  // Creates cycle A->B->C->A
    g.add_edge(c, d, weight=4) |> ignore
    
    let sccs = g.strongly_connected_components()
    println("Number of strongly connected components: \{sccs.length()}")

    // Process each SCC
    for scc in sccs {
        // Process nodes in each strongly connected component
        println("SCC contains \{scc.length()} nodes")
    }
}
```

#### Finding Bridges and Cut Vertices

```moonbit
test {
    // Create a graph to test bridges and articulation points
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=2) |> ignore
    g.add_edge(c, d, weight=3) |> ignore
    g.add_edge(a, c, weight=4) |> ignore  // Creates alternative path
    
    // Find all bridges
    let bridges = g.find_bridges()
    println("Found \{bridges.length()} bridges")

    // Find all cut vertices (articulation points)
    let articulation_points = g.find_articulation_points()
    println("Found \{articulation_points.length()} articulation points")
}
```

### Minimum Spanning Tree

#### Kruskal's Algorithm

```moonbit
test {
    // Create a weighted graph for MST
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=4) |> ignore
    g.add_edge(a, c, weight=2) |> ignore
    g.add_edge(b, c, weight=1) |> ignore
    g.add_edge(c, d, weight=3) |> ignore
    g.add_edge(b, d, weight=5) |> ignore
    
    // Get the MST edge set
    let mst = g.kruskal_mst()

    // Calculate total MST weight
    let total_weight = mst.fold(init=0, fn(acc, edge) { acc + edge.2 })
    println("MST total weight: \{total_weight}")
}
```

#### Prim's Algorithm

```moonbit
test {
    // Create a weighted graph for MST
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=4) |> ignore
    g.add_edge(a, c, weight=2) |> ignore
    g.add_edge(b, c, weight=1) |> ignore
    g.add_edge(c, d, weight=3) |> ignore
    g.add_edge(b, d, weight=5) |> ignore
    
    // Compute MST using Prim's algorithm
    let mst = g.prim_mst()

    // Process MST edges
    for i in 0..<mst.length() {
        let (u, v, weight) = mst[i]
        // Process each edge in the MST
        println("MST edge: \{g.node_weights[u.ix]} -> \{g.node_weights[v.ix]} (weight: \{weight})")
    }
}
```

### Flow Network Algorithms

#### Maximum Flow

```moonbit
test {
    // Create a flow network
    let g : Graph2[String, Int] = Graph2::new(directed=true)
    let s = g.add_node(weight="S")  // Source
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let t = g.add_node(weight="T")  // Sink
    g.add_edge(s, a, weight=10) |> ignore
    g.add_edge(s, b, weight=5) |> ignore
    g.add_edge(a, t, weight=8) |> ignore
    g.add_edge(b, t, weight=7) |> ignore
    g.add_edge(a, b, weight=3) |> ignore
    
    // Calculate maximum flow from source to sink
    let source = s
    let sink = t
    let max_flow = g.max_flow(source, sink)
    println("Maximum flow: \{max_flow}")
}
```

#### Minimum Cost Maximum Flow

```moonbit
test {
    // Create a flow network with costs
    let g : Graph2[String, Int] = Graph2::new(directed=true)
    let s = g.add_node(weight="S")  // Source
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let t = g.add_node(weight="T")  // Sink
    g.add_edge(s, a, weight=10) |> ignore
    g.add_edge(s, b, weight=5) |> ignore
    g.add_edge(a, t, weight=8) |> ignore
    g.add_edge(b, t, weight=7) |> ignore
    g.add_edge(a, b, weight=3) |> ignore
    
    // Define edge costs
    let costs = @hashmap.new()
    costs.set((s, a), 2) |> ignore
    costs.set((s, b), 1) |> ignore
    costs.set((a, t), 3) |> ignore
    costs.set((b, t), 2) |> ignore
    costs.set((a, b), 1) |> ignore

    // Calculate minimum cost maximum flow
    let source = s
    let sink = t
    let (flow, cost) = g.min_cost_max_flow(source, sink, costs)
    println("Flow: \{flow}, Cost: \{cost}")
}
```

#### Dinic Maximum Flow Algorithm

```moonbit
test {
    // Create a flow network for Dinic algorithm
    let g : Graph2[String, Int] = Graph2::new(directed=true)
    let source = g.add_node(weight="source")
    let a = g.add_node(weight="a")
    let b = g.add_node(weight="b")
    let c = g.add_node(weight="c")
    let d = g.add_node(weight="d")
    let sink = g.add_node(weight="sink")
    g.add_edge(source, a, weight=10) |> ignore
    g.add_edge(source, b, weight=10) |> ignore
    g.add_edge(a, c, weight=4) |> ignore
    g.add_edge(a, d, weight=8) |> ignore
    g.add_edge(b, c, weight=9) |> ignore
    g.add_edge(b, d, weight=6) |> ignore
    g.add_edge(c, sink, weight=10) |> ignore
    g.add_edge(d, sink, weight=10) |> ignore
    
    // Calculate maximum flow using Dinic algorithm
    let max_flow = g.dinic_max_flow(source, sink)
    println("Maximum flow (Dinic): \{max_flow}")
}
```

#### Push-Relabel Maximum Flow Algorithm

```moonbit
test {
    // Create a flow network for Push-Relabel algorithm
    let g : Graph2[String, Int] = Graph2::new(directed=true)
    let source = g.add_node(weight="source")
    let a = g.add_node(weight="a")
    let b = g.add_node(weight="b")
    let c = g.add_node(weight="c")
    let d = g.add_node(weight="d")
    let sink = g.add_node(weight="sink")
    g.add_edge(source, a, weight=10) |> ignore
    g.add_edge(source, b, weight=10) |> ignore
    g.add_edge(a, c, weight=4) |> ignore
    g.add_edge(a, d, weight=8) |> ignore
    g.add_edge(b, c, weight=9) |> ignore
    g.add_edge(b, d, weight=6) |> ignore
    g.add_edge(c, sink, weight=10) |> ignore
    g.add_edge(d, sink, weight=10) |> ignore
    
    // Calculate maximum flow using Push-Relabel algorithm
    let max_flow = g.push_relabel_max_flow(source, sink)
    println("Maximum flow (Push-Relabel): \{max_flow}")
}
```

#### Minimum Cost Flow

```moonbit
test {
    // Create a flow network with costs
    let g : Graph2[String, Int] = Graph2::new(directed=true)
    let s = g.add_node(weight="s")
    let a = g.add_node(weight="a")
    let b = g.add_node(weight="b")
    let t = g.add_node(weight="t")
    g.add_edge(s, a, weight=5) |> ignore
    g.add_edge(s, b, weight=4) |> ignore
    g.add_edge(a, b, weight=3) |> ignore
    g.add_edge(a, t, weight=3) |> ignore
    g.add_edge(b, t, weight=8) |> ignore
    
    // Define costs
    let costs = @hashmap.new()
    costs.set((s, a), 2) |> ignore
    costs.set((s, b), 1) |> ignore
    costs.set((a, b), 2) |> ignore
    costs.set((a, t), 3) |> ignore
    costs.set((b, t), 2) |> ignore
    
    // Calculate minimum cost flow with target flow of 7
    let target_flow = Some(7)
    let (flow, cost) = g.min_cost_flow(s, t, target_flow, costs)
    println("Flow: \{flow}, Cost: \{cost}")
}
```

#### Minimum Feasible Flow

```moonbit
test {
    // Create a flow network with lower bounds
    let g : Graph2[String, Int] = Graph2::new(directed=true)
    let s = g.add_node(weight="s")
    let a = g.add_node(weight="a")
    let b = g.add_node(weight="b")
    let t = g.add_node(weight="t")
    g.add_edge(s, a, weight=5) |> ignore
    g.add_edge(a, b, weight=5) |> ignore
    g.add_edge(b, t, weight=5) |> ignore
    g.add_edge(s, b, weight=3) |> ignore
    g.add_edge(a, t, weight=2) |> ignore
    
    // Define lower bounds for edges
    let lower_bounds = Array::new()
    lower_bounds.push((s, a, 1))  // s->a must have at least flow 1
    lower_bounds.push((a, b, 2))  // a->b must have at least flow 2
    lower_bounds.push((b, t, 3))  // b->t must have at least flow 3
    
    // Calculate minimum feasible flow
    let min_flow = g.min_flow(s, t, lower_bounds)
    if min_flow != None {
        println("Minimum feasible flow: \{min_flow.unwrap()}")
    } else {
        println("No feasible flow exists")
    }
}
```

### Graph Feature Analysis

#### Cycle Detection

```moonbit
test {
    // Create a graph with cycles
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=2) |> ignore
    g.add_edge(c, a, weight=3) |> ignore  // Creates cycle A->B->C->A
    g.add_edge(c, d, weight=4) |> ignore
    
    // Find all cycles in the graph
    let cycles = g.find_cycles()
    println("Found \{cycles.length()} cycles")

    // Check each cycle
    for cycle in cycles {
        // Process each cycle
        println("Cycle length: \{cycle.length()}")
    }
}
```

#### Bipartite Graph Detection

```moonbit
test {
    // Create a bipartite graph
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, c, weight=1) |> ignore
    g.add_edge(a, d, weight=2) |> ignore
    g.add_edge(b, c, weight=3) |> ignore
    g.add_edge(b, d, weight=4) |> ignore
    
    // Check if graph is bipartite
    let (is_bipartite, colors) = g.is_bipartite()
    if is_bipartite {
        println("Graph is bipartite")
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
        println("Left part: \{left.length()} nodes, Right part: \{right.length()} nodes")
    } else {
        println("Graph is not bipartite")
    }
}
```

#### Euler Path/Circuit

```moonbit
test {
    // Create a graph for Euler path/circuit testing
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=2) |> ignore
    g.add_edge(c, d, weight=3) |> ignore
    g.add_edge(d, a, weight=4) |> ignore  // Creates Euler circuit
    
    // Check for Euler circuit (all vertices have even degree)
    let has_circuit = g.has_euler_circuit()
    println("Has Euler circuit: \{has_circuit}")

    // Check for Euler path (at most two vertices with odd degree)
    let has_path = g.has_euler_path()
    println("Has Euler path: \{has_path}")
}
```

### Graph Metrics

#### Degree Centrality

```moonbit
test {
    // Create a graph for centrality analysis
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(a, c, weight=2) |> ignore
    g.add_edge(a, d, weight=3) |> ignore  // A has highest degree
    g.add_edge(b, c, weight=4) |> ignore
    
    // Calculate degree centrality for each node
    let centrality = g.degree_centrality()
    println("Degree centrality: \{centrality}")

    // Find the node with highest centrality
    let mut max_centrality = 0.0
    let mut max_node = -1
    for i in 0..<centrality.length() {
        if centrality[i].to_double() > max_centrality {
            max_centrality = centrality[i].to_double()
            max_node = i
        }
    }
    println("Node with highest centrality: \{max_node} (\{max_centrality})")
}
```

#### PageRank

```moonbit
test {
    // Create a directed graph for PageRank
    let g : Graph2[String, Int] = Graph2::new(directed=true)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(a, c, weight=2) |> ignore
    g.add_edge(b, c, weight=3) |> ignore
    g.add_edge(c, d, weight=4) |> ignore
    g.add_edge(d, a, weight=5) |> ignore  // Creates cycle
    
    // Calculate PageRank
    let ranks = g.pagerank(0.85, 100, 0.0001)
    println("PageRank values: \{ranks}")

    // Process PageRank values for each node
    for i in 0..<ranks.length() {
        println("PageRank of node \{i}: \{ranks[i]}")
    }
}
```

#### Clustering Coefficient

```moonbit
test {
    // Create a graph for clustering analysis
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(a, c, weight=2) |> ignore
    g.add_edge(b, c, weight=3) |> ignore  // Creates triangle
    g.add_edge(c, d, weight=4) |> ignore
    
    // Calculate the average clustering coefficient of the graph
    let cc = g.clustering_coefficient()
    println("Clustering coefficient: \{cc}")
}
```

### Graph Coloring

```moonbit
test {
    // Create a graph for coloring
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(a, c, weight=2) |> ignore
    g.add_edge(b, c, weight=3) |> ignore
    g.add_edge(c, d, weight=4) |> ignore
    
    // Color the graph using greedy algorithm
    let (color_count, colors) = g.greedy_coloring()
    println("Greedy coloring: \{color_count} colors, colors: \{colors}")

    // Color the graph using DSatur algorithm
    let (color_count2, colors2) = g.dsatur_coloring()
    println("DSatur coloring: \{color_count2} colors, colors: \{colors2}")
}
```

## Advanced Usage

### Compressed Graph Representation

```moonbit
test {
    // Create a graph
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=2) |> ignore
    g.add_edge(c, d, weight=3) |> ignore
    g.add_edge(a, d, weight=4) |> ignore
    
    // Create compressed graph representation to optimize memory usage
    let compressed = g.to_compressed()
    println("Compressed graph created")

    // Use the compressed graph
    let node = a
    let neighbors = compressed.get_neighbors(node)
    println("Neighbors of node A: \{neighbors.length()}")
}
```

### Graph Diameter

```moonbit
test {
    // Create a tree for diameter calculation
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    let e = g.add_node(weight="E")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=2) |> ignore
    g.add_edge(c, d, weight=3) |> ignore
    g.add_edge(d, e, weight=4) |> ignore
    
    // Calculate tree diameter (longest path)
    let (diameter, path) = g.tree_diameter()
    println("Tree diameter: \{diameter}")
    println("Diameter path length: \{path.length()}")
}
```

## Notes

- Nodes and edges can have weights of various types
- After deleting nodes, node indices may change and need careful handling
- Use `|> ignore` to handle function calls that don't need return values

## Advanced Graph Algorithms

### Graph Isomorphism Detection

```moonbit
test {
    // Create two graphs to test isomorphism
    let g1 : Graph2[Int, Unit] = Graph2::new(directed=false)
    let g2 : Graph2[Int, Unit] = Graph2::new(directed=false)
    
    // Build g1: 0--1--2--3
    //          |     |
    //          4-----5
    let a1 = g1.add_node(weight=0)
    let b1 = g1.add_node(weight=1)
    let c1 = g1.add_node(weight=2)
    let d1 = g1.add_node(weight=3)
    let e1 = g1.add_node(weight=4)
    let f1 = g1.add_node(weight=5)
    g1.add_edge(a1, b1) |> ignore
    g1.add_edge(b1, c1) |> ignore
    g1.add_edge(c1, d1) |> ignore
    g1.add_edge(a1, e1) |> ignore
    g1.add_edge(e1, f1) |> ignore
    g1.add_edge(c1, f1) |> ignore
    
    // Build g2: 0--2--3--5
    //          |     |
    //          1-----4
    let a2 = g2.add_node(weight=0)
    let b2 = g2.add_node(weight=1)
    let c2 = g2.add_node(weight=2)
    let d2 = g2.add_node(weight=3)
    let e2 = g2.add_node(weight=4)
    let f2 = g2.add_node(weight=5)
    g2.add_edge(a2, c2) |> ignore
    g2.add_edge(c2, d2) |> ignore
    g2.add_edge(d2, f2) |> ignore
    g2.add_edge(a2, b2) |> ignore
    g2.add_edge(b2, e2) |> ignore
    g2.add_edge(d2, e2) |> ignore
    
    // Test if graphs are isomorphic
    let is_isomorphic = g1.is_isomorphic(g2)
    println("Graphs are isomorphic: \{is_isomorphic}")
}
```

### Maximum Clique Detection

```moonbit
test {
    // Create a graph with multiple cliques
    let g : Graph2[String, Unit] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    let e = g.add_node(weight="E")
    g.add_edge(a, b) |> ignore
    g.add_edge(a, c) |> ignore
    g.add_edge(b, c) |> ignore  // Triangle ABC
    g.add_edge(c, d) |> ignore
    g.add_edge(c, e) |> ignore
    g.add_edge(d, e) |> ignore  // Triangle CDE
    
    // Find all maximum cliques
    let cliques = g.find_maximum_cliques()
    println("Found \{cliques.length()} maximum cliques")
    
    for i in 0..<cliques.length() {
        let clique = cliques[i]
        println("Clique \{i + 1}: \{clique.length()} nodes")
    }
}
```

### Second Minimum Spanning Tree

```moonbit
test {
    // Create a graph for second MST calculation
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=2) |> ignore
    g.add_edge(c, d, weight=3) |> ignore
    g.add_edge(d, a, weight=4) |> ignore
    g.add_edge(a, c, weight=5) |> ignore
    
    // Calculate second minimum spanning tree weight
    let second_mst_weight = g.second_mst()
    if second_mst_weight != None {
        println("Second MST weight: \{second_mst_weight.unwrap()}")
    } else {
        println("No second MST exists")
    }
}
```

### Minimum Bottleneck Path

```moonbit
test {
    // Create a graph for bottleneck path analysis
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    let e = g.add_node(weight="E")
    let f = g.add_node(weight="F")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=5) |> ignore
    g.add_edge(a, d, weight=3) |> ignore
    g.add_edge(b, e, weight=2) |> ignore
    g.add_edge(c, f, weight=4) |> ignore
    g.add_edge(d, e, weight=6) |> ignore
    g.add_edge(e, f, weight=7) |> ignore
    
    // Find minimum bottleneck path from A to F
    let bottleneck = g.min_bottleneck_path(a, f)
    if bottleneck != None {
        println("Minimum bottleneck: \{bottleneck.unwrap()}")
    } else {
        println("No path exists")
    }
}
```

### Euler Circuit and Path Finding

```moonbit
test {
    // Create a graph with Euler circuit
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=2) |> ignore
    g.add_edge(c, d, weight=3) |> ignore
    g.add_edge(d, a, weight=4) |> ignore
    
    // Check for Euler circuit
    let has_circuit = g.has_euler_circuit()
    println("Has Euler circuit: \{has_circuit}")
    
    // Find Euler circuit if it exists
    let circuit = g.find_euler_circuit()
    if circuit.length() > 0 {
        println("Euler circuit found with \{circuit.length()} nodes")
    } else {
        println("No Euler circuit exists")
    }
    
    // Check for Euler path
    let has_path = g.has_euler_path()
    println("Has Euler path: \{has_path}")
    
    // Find Euler path if it exists
    let path = g.find_euler_path()
    if path.length() > 0 {
        println("Euler path found with \{path.length()} nodes")
    } else {
        println("No Euler path exists")
    }
}
```

### Graph Density Analysis

```moonbit
test {
    // Create a complete graph (high density)
    let g1 : Graph2[Unit, Unit] = Graph2::new(directed=false)
    let a1 = g1.add_node()
    let b1 = g1.add_node()
    let c1 = g1.add_node()
    let d1 = g1.add_node()
    g1.add_edge(a1, b1) |> ignore
    g1.add_edge(a1, c1) |> ignore
    g1.add_edge(a1, d1) |> ignore
    g1.add_edge(b1, c1) |> ignore
    g1.add_edge(b1, d1) |> ignore
    g1.add_edge(c1, d1) |> ignore
    
    let density1 = g1.graph_density()
    println("Complete graph density: \{density1}")
    
    // Create a sparse graph (low density)
    let g2 : Graph2[Unit, Unit] = Graph2::new(directed=false)
    let a2 = g2.add_node()
    let b2 = g2.add_node()
    let c2 = g2.add_node()
    let d2 = g2.add_node()
    let e2 = g2.add_node()
    g2.add_edge(a2, b2) |> ignore
    g2.add_edge(b2, c2) |> ignore
    g2.add_edge(c2, d2) |> ignore
    g2.add_edge(d2, e2) |> ignore
    
    let density2 = g2.graph_density()
    println("Sparse graph density: \{density2}")
}
```

### Node Removal and Graph Modification

```moonbit
test {
    // Create a graph for node removal testing
    let g : Graph2[String, Int] = Graph2::new(directed=false)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=2) |> ignore
    g.add_edge(a, c, weight=3) |> ignore
    
    println("Original graph: \{g.node_count()} nodes, \{g.edge_count()} edges")
    
    // Remove a node and its incident edges
    let removed_weight = g.remove_node(b)
    if removed_weight != None {
        println("Removed node with weight: \{removed_weight.unwrap()}")
    }
    
    println("After removal: \{g.node_count()} nodes, \{g.edge_count()} edges")
}
```

## Advanced Examples

### Complex Network Analysis

```moonbit
test {
    // Create a complex network
    let g : Graph2[String, Int] = Graph2::new(directed=true)
    let a = g.add_node(weight="A")
    let b = g.add_node(weight="B")
    let c = g.add_node(weight="C")
    let d = g.add_node(weight="D")
    let e = g.add_node(weight="E")
    let f = g.add_node(weight="F")
    
    // Add nodes and edges to create a complex network
    g.add_edge(a, b, weight=1) |> ignore
    g.add_edge(b, c, weight=2) |> ignore
    g.add_edge(c, a, weight=3) |> ignore  // Cycle 1
    g.add_edge(d, e, weight=4) |> ignore
    g.add_edge(e, f, weight=5) |> ignore
    g.add_edge(f, d, weight=6) |> ignore  // Cycle 2
    g.add_edge(c, d, weight=7) |> ignore  // Connection between cycles

    // Compute strongly connected components
    let sccs = g.strongly_connected_components()
    println("Found \{sccs.length()} strongly connected components")

    // Calculate PageRank for each SCC
    for scc in sccs {
        println("SCC contains \{scc.length()} nodes")
        // Note: subgraph method may not be available in this version
        // Analyze PageRank results for each subgraph
    }
}
```