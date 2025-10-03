///| tests for ds.mbt

test "Graph2 add_node_and_edge" {
  let g: Graph2[Unit, Unit] = Graph2::new(directed=false)
  let a = g.add_node()
  let b = g.add_node()
  g.add_edge(a, b) |> ignore
  inspect(g.node_count().to_string(), content="2")
}

// connected components
test "connected_components" {
  let g: Graph2[Unit, Unit] = Graph2::new(directed=false)
  let a = g.add_node()
  let b = g.add_node()
  let c = g.add_node()
  let d = g.add_node()
  g.add_edge(a, b) |> ignore
  g.add_edge(c, d) |> ignore
  let comps = g.connected_components()
  inspect(comps.length().to_string(), content="2")
}

// strongly connected components (Kosaraju)
test "strongly_connected_components" {
  let g: Graph2[Unit, Unit] = Graph2::new(directed=true)
  // build cycle 0->1->2->0 and a tail 3->4
  let n0 = g.add_node()
  let n1 = g.add_node()
  let n2 = g.add_node()
  let n3 = g.add_node()
  let n4 = g.add_node()
  g.add_edge(n0, n1) |> ignore
  g.add_edge(n1, n2) |> ignore
  g.add_edge(n2, n0) |> ignore
  g.add_edge(n3, n4) |> ignore
  let sccs = g.strongly_connected_components()
  // expect one SCC of size 3 and two singletons
  let sizes = sccs.map(fn(c) { c.length() })
  // bubble sort sizes for assertion
  let ss = sizes
  let len_ss = ss.length()
  let mut i = 0
  while i < len_ss {
    let mut j = i + 1
    while j < len_ss {
      if ss[i] > ss[j] {
        ss.swap(i, j)
      }
      j = j + 1
    }
    i = i + 1
  }
  inspect(ss[0].to_string(), content="1")
  inspect(ss[1].to_string(), content="1")
  inspect(ss[2].to_string(), content="3")
}

// dijkstra int
test "dijkstra_int" {
  let g: Graph2[Int, Int] = Graph2::new(directed=false)
  let a = g.add_node(weight=0)
  let b = g.add_node(weight=0)
  let c = g.add_node(weight=0)
  g.add_edge(a, b, weight=1) |> ignore
  g.add_edge(b, c, weight=2) |> ignore
  let dist = g.dijkstra_int(a)
  inspect(dist[b.ix].to_string(), content="1")
  inspect(dist[c.ix].to_string(), content="3")
}

// kruskal mst
test "kruskal_mst" {
  let g: Graph2[Unit, Int] = Graph2::new(directed=false)
  let a = g.add_node()
  let b = g.add_node()
  let c = g.add_node()
  g.add_edge(a, b, weight=1) |> ignore
  g.add_edge(b, c, weight=2) |> ignore
  g.add_edge(a, c, weight=5) |> ignore
  let mst = g.kruskal_mst()
  // mst should have 2 edges with total weight 3
  let mut total = 0
  for e in mst {
    total = total + e.2
  }
  inspect(total.to_string(), content="3")
}
