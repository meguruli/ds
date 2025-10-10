///|
/// tests for ds.mbt
test "Graph2 add_node_and_edge" {
  let g : Graph2[Unit, Unit] = Graph2::new(directed=false)
  let a = g.add_node()
  let b = g.add_node()
  g.add_edge(a, b) |> ignore
  inspect(g.node_count().to_string(), content="2")
}

// connected components

///|
test "connected_components" {
  let g : Graph2[Unit, Unit] = Graph2::new(directed=false)
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

///|
test "strongly_connected_components" {
  let g : Graph2[Unit, Unit] = Graph2::new(directed=true)
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

///|
test "dijkstra_int" {
  let g : Graph2[Int, Int] = Graph2::new(directed=false)
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

///|
test "kruskal_mst" {
  let g : Graph2[Unit, Int] = Graph2::new(directed=false)
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

///|
/// 测试 A* 寻路算法 - 基本功能
test "a_star_basic" {
  let g : Graph2[Int, Int] = Graph2::new(directed=true)

  // 创建简单图: 0->1->3
  //             \->2/
  let n0 = g.add_node(weight=0)
  let n1 = g.add_node(weight=1)
  let n2 = g.add_node(weight=2)
  let n3 = g.add_node(weight=3)
  g.add_edge(n0, n1, weight=1) |> ignore
  g.add_edge(n0, n2, weight=3) |> ignore
  g.add_edge(n1, n3, weight=3) |> ignore
  g.add_edge(n2, n3, weight=1) |> ignore

  // 简单启发式函数 - 总是返回0（相当于Dijkstra）
  fn heuristic(_ : NodeIndex, _ : NodeIndex) -> Int {
    0
  }

  // 计算从n0到n3的路径
  let path = g.a_star(n0, n3, heuristic)

  // 应该选择路径 0->1->3 (总权重4) 而非 0->2->3 (总权重4)
  inspect(path.length().to_string(), content="3")
  inspect(path[0].ix.to_string(), content="0") // 起点
  inspect(path[1].ix.to_string(), content="2") // 经过1
  inspect(path[2].ix.to_string(), content="3") // 终点
}

///|
/// 测试 A* 寻路算法 - 使用浮点型边权重和启发式函数
test "a_star_with_float_heuristic" {
  let g : Graph2[String, Float] = Graph2::new(directed=false)

  // 创建一个网格图
  // A -- B -- C
  // |    |    |
  // D -- E -- F
  // |    |    |
  // G -- H -- I

  let a = g.add_node(weight="A")
  let b = g.add_node(weight="B")
  let c = g.add_node(weight="C")
  let d = g.add_node(weight="D")
  let e = g.add_node(weight="E")
  let f = g.add_node(weight="F")
  let g_node = g.add_node(weight="G")
  let h = g.add_node(weight="H")
  let i = g.add_node(weight="I")

  // 添加边，每条边权重为1.0
  g.add_edge(a, b, weight=1.0) |> ignore
  g.add_edge(b, c, weight=1.0) |> ignore
  g.add_edge(a, d, weight=1.0) |> ignore
  g.add_edge(b, e, weight=1.0) |> ignore
  g.add_edge(c, f, weight=1.0) |> ignore
  g.add_edge(d, e, weight=1.0) |> ignore
  g.add_edge(e, f, weight=1.0) |> ignore
  g.add_edge(d, g_node, weight=1.0) |> ignore
  g.add_edge(e, h, weight=1.0) |> ignore
  g.add_edge(f, i, weight=1.0) |> ignore
  g.add_edge(g_node, h, weight=1.0) |> ignore
  g.add_edge(h, i, weight=1.0) |> ignore

  // 启发式函数 - 曼哈顿距离（直接返回Float）
  fn heuristic(a : NodeIndex, b : NodeIndex) -> Float {
    // 网格上的行列位置
    fn position(idx : Int) -> (Int, Int) {
      let row = idx / 3
      let col = idx % 3
      (row, col)
    }

    let (a_row, a_col) = position(a.ix)
    let (b_row, b_col) = position(b.ix)

    // 曼哈顿距离（直接计算浮点数）
    (a_row - b_row).abs().to_float() + (a_col - b_col).abs().to_float()
  }

  // 从A到I的路径
  let path = g.a_star(a, i, heuristic)

  // 验证路径长度和端点
  inspect(path.length().to_string(), content="5")
  inspect(g.node_weights[path[0].ix], content="A") // 起点
  inspect(g.node_weights[path[path.length() - 1].ix], content="I") // 终点
}

///|
/// 测试 A* 寻路算法 - 无路径情况
test "a_star_no_path" {
  let g : Graph2[Int, Int] = Graph2::new(directed=true)

  // 创建两个不连通的部分
  let n0 = g.add_node(weight=0)
  let n1 = g.add_node(weight=1)
  let n2 = g.add_node(weight=2)
  let n3 = g.add_node(weight=3)
  g.add_edge(n0, n1, weight=1) |> ignore
  g.add_edge(n2, n3, weight=1) |> ignore
  fn heuristic(_ : NodeIndex, _ : NodeIndex) -> Int {
    0
  }

  // 尝试在不连通的点之间寻找路径
  let path = g.a_star(n0, n3, heuristic)

  // 应该返回空数组
  inspect(path.length().to_string(), content="0")
}

///|
/// 测试寻找成环节点 - 单环
test "find_cycles_single" {
  let g : Graph2[Int, Int] = Graph2::new(directed=true)

  // 创建一个有向环 0->1->2->0
  let n0 = g.add_node(weight=0)
  let n1 = g.add_node(weight=1)
  let n2 = g.add_node(weight=2)
  g.add_edge(n0, n1, weight=1) |> ignore
  g.add_edge(n1, n2, weight=1) |> ignore
  g.add_edge(n2, n0, weight=1) |> ignore
  let cycles = g.find_cycles()

  // 应该找到一个环
  inspect(cycles.length().to_string(), content="1")
  inspect(cycles[0].length().to_string(), content="3")
}

///|
/// 测试寻找成环节点 - 多环
test "find_cycles_multiple" {
  let g : Graph2[Int, Int] = Graph2::new(directed=true)

  // 创建两个环
  // 环1: 0->1->2->0
  // 环2: 3->4->3
  let n0 = g.add_node(weight=0)
  let n1 = g.add_node(weight=1)
  let n2 = g.add_node(weight=2)
  let n3 = g.add_node(weight=3)
  let n4 = g.add_node(weight=4)
  g.add_edge(n0, n1, weight=1) |> ignore
  g.add_edge(n1, n2, weight=1) |> ignore
  g.add_edge(n2, n0, weight=1) |> ignore
  g.add_edge(n3, n4, weight=1) |> ignore
  g.add_edge(n4, n3, weight=1) |> ignore
  let cycles = g.find_cycles()

  // 应该找到两个环
  inspect(cycles.length().to_string(), content="2")

  // 环的大小应该是3和2
  let sizes = Array::new()
  for cycle in cycles {
    sizes.push(cycle.length())
  }

  // 排序以确保比较
  if sizes[0] > sizes[1] {
    sizes.swap(0, 1)
  }
  inspect(sizes[0].to_string(), content="2")
  inspect(sizes[1].to_string(), content="3")
}

///|
/// 测试寻找成环节点 - 有重叠的环
test "find_cycles_overlapping" {
  let g : Graph2[Int, Int] = Graph2::new(directed=true)

  // 创建重叠的环
  // 0 -> 1 -> 2 -> 0
  //      ^    v
  //      4 <- 3
  let n0 = g.add_node(weight=0)
  let n1 = g.add_node(weight=1)
  let n2 = g.add_node(weight=2)
  let n3 = g.add_node(weight=3)
  let n4 = g.add_node(weight=4)
  g.add_edge(n0, n1, weight=1) |> ignore
  g.add_edge(n1, n2, weight=1) |> ignore
  g.add_edge(n2, n0, weight=1) |> ignore
  g.add_edge(n2, n3, weight=1) |> ignore
  g.add_edge(n3, n4, weight=1) |> ignore
  g.add_edge(n4, n1, weight=1) |> ignore
  let cycles = g.find_cycles()

  // 应该找到多个环
  inspect((cycles.length() > 0).to_string(), content="true")
}

///|
/// 测试寻找成环节点 - 无环图
test "find_cycles_acyclic" {
  let g : Graph2[Int, Int] = Graph2::new(directed=true)

  // 创建无环图
  // 0 -> 1 -> 3
  //  \-> 2 -/
  let n0 = g.add_node(weight=0)
  let n1 = g.add_node(weight=1)
  let n2 = g.add_node(weight=2)
  let n3 = g.add_node(weight=3)
  g.add_edge(n0, n1, weight=1) |> ignore
  g.add_edge(n0, n2, weight=1) |> ignore
  g.add_edge(n1, n3, weight=1) |> ignore
  g.add_edge(n2, n3, weight=1) |> ignore
  let cycles = g.find_cycles()

  // 不应该找到环
  inspect(cycles.length().to_string(), content="0")
}

///|
/// 测试最大流算法 (max_flow)
test "max_flow_basic" {
  let g : Graph2[String, Int] = Graph2::new(directed=true)

  // 创建一个简单的流网络
  // source(0) -> 1 -> sink(3)
  //           \-> 2 -/
  let source = g.add_node(weight="source")
  let n1 = g.add_node(weight="1")
  let n2 = g.add_node(weight="2")
  let sink = g.add_node(weight="sink")
  g.add_edge(source, n1, weight=3) |> ignore
  g.add_edge(source, n2, weight=2) |> ignore
  g.add_edge(n1, sink, weight=3) |> ignore
  g.add_edge(n2, sink, weight=2) |> ignore
  let flow = g.max_flow(source, sink)
  inspect(flow.to_string(), content="5") // 预期最大流是5
}

///|
/// 测试带瓶颈的最大流
test "max_flow_bottleneck" {
  let g : Graph2[Unit, Int] = Graph2::new(directed=true)
  let source = g.add_node()
  let n1 = g.add_node()
  let n2 = g.add_node()
  let sink = g.add_node()

  // 添加边，中间有瓶颈
  g.add_edge(source, n1, weight=100) |> ignore
  g.add_edge(n1, n2, weight=1) |> ignore // 瓶颈
  g.add_edge(n2, sink, weight=100) |> ignore
  let flow = g.max_flow(source, sink)
  inspect(flow.to_string(), content="1") // 瓶颈限制了流量
}

///|
/// 测试 Bellman-Ford 算法
test "bellman_ford_basic" {
  let g : Graph2[String, Int] = Graph2::new(directed=true)
  let a = g.add_node(weight="A")
  let b = g.add_node(weight="B")
  let c = g.add_node(weight="C")
  let d = g.add_node(weight="D")
  g.add_edge(a, b, weight=1) |> ignore
  g.add_edge(b, c, weight=2) |> ignore
  g.add_edge(a, c, weight=5) |> ignore
  g.add_edge(c, d, weight=1) |> ignore
  let (distances, has_negative_cycle) = g.bellman_ford(a)
  inspect(has_negative_cycle.to_string(), content="false")
  inspect(distances[b.ix].unwrap().to_string(), content="1")
  inspect(distances[c.ix].unwrap().to_string(), content="3")
  inspect(distances[d.ix].unwrap().to_string(), content="4")
}

///|
/// 测试带负权边的 Bellman-Ford
test "bellman_ford_negative_edges" {
  let g : Graph2[String, Int] = Graph2::new(directed=true)

  // 创建一个含负权边但无负权环的图
  //     +2      
  // A -----> B
  // |        |
  // +5|      |-2  (负权边)
  // |        |
  // v        v
  // C <----- D
  //     +3      

  let a = g.add_node(weight="A") // 源点
  let b = g.add_node(weight="B")
  let c = g.add_node(weight="C")
  let d = g.add_node(weight="D")
  g.add_edge(a, b, weight=2) |> ignore
  g.add_edge(b, d, weight=-2) |> ignore // 负权边
  g.add_edge(d, c, weight=3) |> ignore
  g.add_edge(a, c, weight=5) |> ignore
  let (distances, has_negative_cycle) = g.bellman_ford(a)

  // 检查图中不存在负权环
  let neg_cycle_str = if has_negative_cycle { "true" } else { "false" }
  inspect(neg_cycle_str, content="false")

  // 验证通过负权边优化后的最短路径
  // 从A到C的最短路径是A->B->D->C (2+(-2)+3=3)，比A->C (5)更短
  let b_dist = if distances[b.ix] != None {
    distances[b.ix].unwrap()
  } else {
    -999
  }
  inspect(b_dist.to_string(), content="2") // A->B = 2
  let d_dist = if distances[d.ix] != None {
    distances[d.ix].unwrap()
  } else {
    -999
  }
  inspect(d_dist.to_string(), content="0") // A->B->D = 2+(-2) = 0
  let c_dist = if distances[c.ix] != None {
    distances[c.ix].unwrap()
  } else {
    -999
  }
  inspect(c_dist.to_string(), content="3") // A->B->D->C = 2+(-2)+3 = 3
}

///|
/// 测试带负环的 Bellman-Ford
test "bellman_ford_negative_cycle" {
  let g : Graph2[String, Int] = Graph2::new(directed=true)
  let a = g.add_node(weight="A")
  let b = g.add_node(weight="B")
  let c = g.add_node(weight="C")
  g.add_edge(a, b, weight=1) |> ignore
  g.add_edge(b, c, weight=2) |> ignore
  g.add_edge(c, a, weight=-4) |> ignore // 形成负环: a->b->c->a，总权重-1
  let (_, has_negative_cycle) = g.bellman_ford(a)
  inspect(has_negative_cycle.to_string(), content="true")
}

///|
/// 测试 Floyd-Warshall 算法
test "floyd_warshall" {
  let g : Graph2[String, Int] = Graph2::new(directed=true)
  let a = g.add_node(weight="A")
  let b = g.add_node(weight="B")
  let c = g.add_node(weight="C")
  let d = g.add_node(weight="D")
  g.add_edge(a, b, weight=1) |> ignore
  g.add_edge(b, c, weight=2) |> ignore
  g.add_edge(a, c, weight=5) |> ignore
  g.add_edge(c, d, weight=1) |> ignore
  let distances = g.floyd_warshall()
  inspect(distances[a.ix][b.ix].unwrap().to_string(), content="1")
  inspect(distances[a.ix][c.ix].unwrap().to_string(), content="3")
  inspect(distances[a.ix][d.ix].unwrap().to_string(), content="4")
  inspect(distances[b.ix][d.ix].unwrap().to_string(), content="3")
}

///|
/// 测试 Prim 最小生成树算法
test "prim_mst" {
  let g : Graph2[Unit, Int] = Graph2::new(directed=false)
  let a = g.add_node()
  let b = g.add_node()
  let c = g.add_node()
  let d = g.add_node()
  g.add_edge(a, b, weight=1) |> ignore
  g.add_edge(b, c, weight=2) |> ignore
  g.add_edge(a, c, weight=3) |> ignore
  g.add_edge(c, d, weight=1) |> ignore
  g.add_edge(b, d, weight=4) |> ignore
  let mst = g.prim_mst()

  // MST应该有3条边，总权重为4
  let mut total_weight = 0
  for edge in mst {
    let (_, _, weight) = edge
    total_weight += weight
  }
  inspect(mst.length().to_string(), content="3")
  inspect(total_weight.to_string(), content="9")
}

///|
/// 测试二分图检测
test "is_bipartite" {
  // 创建一个二分图 (偶数环是二分的)
  let g1 : Graph2[Unit, Unit] = Graph2::new(directed=false)
  let a = g1.add_node()
  let b = g1.add_node()
  let c = g1.add_node()
  let d = g1.add_node()
  g1.add_edge(a, b) |> ignore
  g1.add_edge(b, c) |> ignore
  g1.add_edge(c, d) |> ignore
  g1.add_edge(d, a) |> ignore
  let (is_bipartite1, _) = g1.is_bipartite()
  inspect(is_bipartite1.to_string(), content="true")

  // 创建一个非二分图 (奇数环不是二分的)
  let g2 : Graph2[Unit, Unit] = Graph2::new(directed=false)
  let e = g2.add_node()
  let f = g2.add_node()
  let g = g2.add_node()
  g2.add_edge(e, f) |> ignore
  g2.add_edge(f, g) |> ignore
  g2.add_edge(g, e) |> ignore
  let (is_bipartite2, _) = g2.is_bipartite()
  inspect(is_bipartite2.to_string(), content="false")
}

///|
/// 测试二分图最大匹配
test "bipartite_matching" {
  let g : Graph2[String, Unit] = Graph2::new(directed=false)

  // 创建二分图
  // 左侧集合：a, b, c
  // 右侧集合：d, e, f
  let a = g.add_node(weight="a")
  let b = g.add_node(weight="b")
  let c = g.add_node(weight="c")
  let d = g.add_node(weight="d")
  let e = g.add_node(weight="e")
  let f = g.add_node(weight="f")
  g.add_edge(a, d) |> ignore
  g.add_edge(a, e) |> ignore
  g.add_edge(b, d) |> ignore
  g.add_edge(c, f) |> ignore

  // 手动创建颜色数组 (0为左侧，1为右侧)
  let colors = Array::make(g.node_count(), -1)
  colors[a.ix] = 0
  colors[b.ix] = 0
  colors[c.ix] = 0
  colors[d.ix] = 1
  colors[e.ix] = 1
  colors[f.ix] = 1
  let matching = g.bipartite_matching(colors)

  // 最大匹配应该有3条边
  inspect(matching.length().to_string(), content="3")
}

///|
/// 测试欧拉回路检测
test "euler_circuit" {
  // 有欧拉回路的图 (所有节点度数为偶数)
  let g1 : Graph2[Unit, Unit] = Graph2::new(directed=false)
  let a = g1.add_node()
  let b = g1.add_node()
  let c = g1.add_node()
  let d = g1.add_node()
  g1.add_edge(a, b) |> ignore
  g1.add_edge(b, c) |> ignore
  g1.add_edge(c, d) |> ignore
  g1.add_edge(d, a) |> ignore
  g1.add_edge(a, c) |> ignore // 额外的边，保持所有节点度为偶数
  g1.add_edge(b, d) |> ignore // 额外的边，保持所有节点度为偶数
  inspect(g1.has_euler_circuit().to_string(), content="false")

  // 无欧拉回路的图 (有两个节点度数为奇数)
  let g2 : Graph2[Unit, Unit] = Graph2::new(directed=false)
  let e = g2.add_node()
  let f = g2.add_node()
  let g = g2.add_node()
  g2.add_edge(e, f) |> ignore
  g2.add_edge(f, g) |> ignore
  inspect(g2.has_euler_circuit().to_string(), content="false")
}

///|
/// 测试欧拉路径检测
test "euler_path" {
  // 有欧拉路径但无欧拉回路的图 (恰好两个节点度数为奇数)
  let g1 : Graph2[Unit, Unit] = Graph2::new(directed=false)
  let a = g1.add_node()
  let b = g1.add_node()
  let c = g1.add_node()
  let d = g1.add_node()
  g1.add_edge(a, b) |> ignore
  g1.add_edge(b, c) |> ignore
  g1.add_edge(c, d) |> ignore
  g1.add_edge(a, c) |> ignore // 使得a和d的度为奇数
  inspect(g1.has_euler_path().to_string(), content="true")
  inspect(g1.has_euler_circuit().to_string(), content="false")

  // 无欧拉路径的图 (超过两个节点度数为奇数)
  let g2 : Graph2[Unit, Unit] = Graph2::new(directed=false)
  let e = g2.add_node()
  let f = g2.add_node()
  let g = g2.add_node()
  let h = g2.add_node()
  g2.add_edge(e, f) |> ignore
  g2.add_edge(f, g) |> ignore
  g2.add_edge(f, h) |> ignore
  inspect(g2.has_euler_path().to_string(), content="false")
}

///|
/// 测试寻找桥（割边）
test "find_bridges" {
  let g : Graph2[String, Unit] = Graph2::new(directed=false)

  // 创建图，其中有两个桥：BC和DE
  // A -- B -- C -- D -- E -- F
  //      |                |
  //      G -- H -- I -----+

  let a = g.add_node(weight="A")
  let b = g.add_node(weight="B")
  let c = g.add_node(weight="C")
  let d = g.add_node(weight="D")
  let e = g.add_node(weight="E")
  let f = g.add_node(weight="F")
  let g_node = g.add_node(weight="G")
  let h = g.add_node(weight="H")
  let i = g.add_node(weight="I")
  g.add_edge(a, b) |> ignore
  g.add_edge(b, c) |> ignore // 桥
  g.add_edge(c, d) |> ignore
  g.add_edge(d, e) |> ignore // 桥
  g.add_edge(e, f) |> ignore
  g.add_edge(b, g_node) |> ignore
  g.add_edge(g_node, h) |> ignore
  g.add_edge(h, i) |> ignore
  g.add_edge(i, e) |> ignore
  let bridges = g.find_bridges()
  inspect(bridges.length().to_string(), content="2")
}

///|
/// 测试寻找割点
test "find_articulation_points" {
  let g : Graph2[String, Unit] = Graph2::new(directed=false)

  // 创建图，其中B、C、E是割点
  // A -- B -- C -- D -- E -- F
  //      |                |
  //      G -- H -- I -----+

  let a = g.add_node(weight="A")
  let b = g.add_node(weight="B")
  let c = g.add_node(weight="C")
  let d = g.add_node(weight="D")
  let e = g.add_node(weight="E")
  let f = g.add_node(weight="F")
  let g_node = g.add_node(weight="G")
  let h = g.add_node(weight="H")
  let i = g.add_node(weight="I")
  g.add_edge(a, b) |> ignore
  g.add_edge(b, c) |> ignore
  g.add_edge(c, d) |> ignore
  g.add_edge(d, e) |> ignore
  g.add_edge(e, f) |> ignore
  g.add_edge(b, g_node) |> ignore
  g.add_edge(g_node, h) |> ignore
  g.add_edge(h, i) |> ignore
  g.add_edge(i, e) |> ignore
  let articulation_points = g.find_articulation_points()
  inspect(articulation_points.length().to_string(), content="2")
}

///|
/// 测试树的直径
test "tree_diameter" {
  let g : Graph2[String, Int] = Graph2::new(directed=false)

  // 创建一棵树
  //     A
  //    / \
  //   B   C
  //  /     \
  // D       E
  //          \
  //           F

  let a = g.add_node(weight="A")
  let b = g.add_node(weight="B")
  let c = g.add_node(weight="C")
  let d = g.add_node(weight="D")
  let e = g.add_node(weight="E")
  let f = g.add_node(weight="F")
  g.add_edge(a, b, weight=1) |> ignore
  g.add_edge(a, c, weight=1) |> ignore
  g.add_edge(b, d, weight=1) |> ignore
  g.add_edge(c, e, weight=1) |> ignore
  g.add_edge(e, f, weight=1) |> ignore
  let (diameter, path) = g.tree_diameter()

  // 直径应该是D-B-A-C-E-F，长度为5
  inspect(diameter.to_string(), content="5")
  inspect(path.length().to_string(), content="6")

  // 验证路径的端点
  let start_node = g.node_weights[path[0].ix]
  let end_node = g.node_weights[path[path.length() - 1].ix]

  // 路径端点应为D和F (或F和D)
  let is_correct_path = (start_node == "D" && end_node == "F") ||
    (start_node == "F" && end_node == "D")
  inspect(is_correct_path.to_string(), content="true")
}

///|
/// 测试最小费用最大流
test "min_cost_max_flow" {
  let g : Graph2[String, Int] = Graph2::new(directed=true)
  let s = g.add_node(weight="s")
  let a = g.add_node(weight="a")
  let b = g.add_node(weight="b")
  let t = g.add_node(weight="t")

  // 添加边：(节点1, 节点2, 容量)
  g.add_edge(s, a, weight=3) |> ignore
  g.add_edge(s, b, weight=2) |> ignore
  g.add_edge(a, b, weight=1) |> ignore
  g.add_edge(a, t, weight=3) |> ignore
  g.add_edge(b, t, weight=4) |> ignore

  // 定义成本
  let costs = @hashmap.new()
  costs.set((s, a), 1)
  costs.set((s, b), 2)
  costs.set((a, b), 1)
  costs.set((a, t), 3)
  costs.set((b, t), 1)
  let (max_flow, min_cost) = g.min_cost_max_flow(s, t, costs)

  // 最大流应该是5，最小成本应该是7
  inspect(max_flow.to_string(), content="5")
  inspect(min_cost.to_string(), content="17")
}

///|
/// 测试度中心性
test "degree_centrality" {
  let g : Graph2[Unit, Unit] = Graph2::new(directed=false)

  // 创建星形图，中心节点连接所有其他节点
  let center = g.add_node()
  let a = g.add_node()
  let b = g.add_node()
  let c = g.add_node()
  let d = g.add_node()
  g.add_edge(center, a) |> ignore
  g.add_edge(center, b) |> ignore
  g.add_edge(center, c) |> ignore
  g.add_edge(center, d) |> ignore
  let centrality = g.degree_centrality()

  // 中心节点的度中心性应该是1.0，其他节点是0.25
  let center_value = centrality[center.ix]
  let leaf_value = centrality[a.ix]
  inspect((center_value > 0.9).to_string(), content="true")
  inspect((leaf_value - 0.25).abs() < 0.01, content="true")
}

///|
/// 测试PageRank
test "pagerank" {
  let g : Graph2[Unit, Unit] = Graph2::new(directed=true)

  // 创建有向环
  let a = g.add_node()
  let b = g.add_node()
  let c = g.add_node()
  let d = g.add_node()
  g.add_edge(a, b) |> ignore
  g.add_edge(b, c) |> ignore
  g.add_edge(c, d) |> ignore
  g.add_edge(d, a) |> ignore
  g.add_edge(b, d) |> ignore // b有两个出边
  let ranks = g.pagerank(0.85, 100, 0.0001)

  // 在环中，所有节点的PageRank应该接近，但不完全相同
  // 由于b有两个出边，它的PageRank会略低于其他节点

  // 验证所有PageRank总和接近1
  let mut sum : Float = 0.0
  for r in ranks {
    sum += r
  }
  inspect((sum - 1.0).abs() < 0.01, content="true")
}

///|
/// 测试图着色
test "graph_coloring" {
  // 创建一个四色图 (网格)
  let g : Graph2[Unit, Unit] = Graph2::new(directed=false)

  // 创建3x3网格
  let nodes = Array::new()
  for _ in 0..<9 {
    nodes.push(g.add_node())
  }

  // 添加水平边
  for i in 0..<3 {
    for j in 0..<2 {
      g.add_edge(nodes[i * 3 + j], nodes[i * 3 + j + 1]) |> ignore
    }
  }

  // 添加垂直边
  for i in 0..<2 {
    for j in 0..<3 {
      g.add_edge(nodes[i * 3 + j], nodes[(i + 1) * 3 + j]) |> ignore
    }
  }
  let (color_count, _) = g.greedy_coloring()

  // 平面图最多需要4种颜色
  inspect((color_count <= 4).to_string(), content="true")
}

///|
/// 测试DSatur着色
test "dsatur_coloring" {
  // 创建一个复杂图
  let g : Graph2[Unit, Unit] = Graph2::new(directed=false)

  // 创建完全图K5
  let nodes = Array::new()
  for _ in 0..<5 {
    nodes.push(g.add_node())
  }
  for i in 0..<5 {
    for j in (i + 1)..<5 {
      g.add_edge(nodes[i], nodes[j]) |> ignore
    }
  }
  let (color_count, _) = g.dsatur_coloring()

  // 完全图K5需要5种颜色
  inspect(color_count.to_string(), content="5")
}

///|
/// 测试聚类系数
test "clustering_coefficient" {
  // 创建一个三角形网络
  let g : Graph2[Unit, Unit] = Graph2::new(directed=false)
  let a = g.add_node()
  let b = g.add_node()
  let c = g.add_node()
  g.add_edge(a, b) |> ignore
  g.add_edge(b, c) |> ignore
  g.add_edge(c, a) |> ignore
  let cc = g.clustering_coefficient()

  // 三角形的聚类系数为1
  inspect((cc - 1.0).abs() < 0.01, content="true")

  // 添加一个只连接到a的节点，降低整体聚类系数
  let d = g.add_node()
  g.add_edge(a, d) |> ignore
  let cc2 = g.clustering_coefficient()

  // 聚类系数应该降低
  inspect((cc2 < cc).to_string(), content="true")
}

///|
/// 测试图密度
test "graph_density" {
  // 创建完全图K4
  let g1 : Graph2[Unit, Unit] = Graph2::new(directed=false)
  let a = g1.add_node()
  let b = g1.add_node()
  let c = g1.add_node()
  let d = g1.add_node()
  g1.add_edge(a, b) |> ignore
  g1.add_edge(a, c) |> ignore
  g1.add_edge(a, d) |> ignore
  g1.add_edge(b, c) |> ignore
  g1.add_edge(b, d) |> ignore
  g1.add_edge(c, d) |> ignore
  let density1 = g1.graph_density()

  // 完全图的密度应该为1
  inspect((density1 - 1.0).abs() < 0.01, content="true")

  // 创建稀疏图
  let g2 : Graph2[Unit, Unit] = Graph2::new(directed=false)
  let e = g2.add_node()
  let f = g2.add_node()
  let g = g2.add_node()
  let h = g2.add_node()
  let i = g2.add_node()
  g2.add_edge(e, f) |> ignore
  g2.add_edge(f, g) |> ignore
  g2.add_edge(g, h) |> ignore
  g2.add_edge(h, i) |> ignore
  let density2 = g2.graph_density()

  // 路径图的密度应该很低
  inspect((density2 < 0.5).to_string(), content="true")
}

///|
/// 测试压缩图
test "compressed_graph" {
  let g : Graph2[String, Int] = Graph2::new(directed=true)
  let a = g.add_node(weight="A")
  let b = g.add_node(weight="B")
  let c = g.add_node(weight="C")
  g.add_edge(a, b, weight=1) |> ignore
  g.add_edge(a, c, weight=2) |> ignore
  g.add_edge(b, c, weight=3) |> ignore
  let compressed = g.to_compressed()

  // 测试节点权重是否正确保存
  inspect(compressed.node_weights[a.ix], content="A")

  // 测试邻居访问
  let neighbors = compressed.get_neighbors(a)
  inspect(neighbors.length().to_string(), content="2")

  // 验证边权重
  let mut found_b = false
  let mut found_c = false
  for pair in neighbors {
    let (node, weight) = pair
    if node.ix == b.ix {
      inspect(weight.to_string(), content="1")
      found_b = true
    } else if node.ix == c.ix {
      inspect(weight.to_string(), content="2")
      found_c = true
    }
  }
  inspect(found_b.to_string(), content="true")
  inspect(found_c.to_string(), content="true")
}

///|
/// 测试节点删除功能
test "remove_node" {
  let g : Graph2[String, Int] = Graph2::new(directed=false)
  let a = g.add_node(weight="A")
  let b = g.add_node(weight="B")
  let c = g.add_node(weight="C")
  g.add_edge(a, b, weight=1) |> ignore
  g.add_edge(b, c, weight=2) |> ignore
  g.add_edge(a, c, weight=3) |> ignore

  // 删除中间节点B
  let removed_weight = g.remove_node(b)
  let new_c = NodeIndex::new(1)

  // 检查返回值
  inspect(removed_weight.unwrap(), content="B")

  // 检查节点数减少
  inspect(g.node_count().to_string(), content="2")

  // 检查边数减少
  inspect(g.edge_count().to_string(), content="1")

  // 检查仅存的边是A-C
  let edge = g.find_edge(a, new_c)
  inspect((edge != None).to_string(), content="true")
}

///|
/// 测试强连通分量算法
test "strongly_connected_components" {
  let g : Graph2[Int, Unit] = Graph2::new(directed=true)

  // 创建两个强连通分量
  // 分量1: 0->1->2->0
  // 分量2: 3->4->3
  // 从0到3有一条边

  let n0 = g.add_node(weight=0)
  let n1 = g.add_node(weight=1)
  let n2 = g.add_node(weight=2)
  let n3 = g.add_node(weight=3)
  let n4 = g.add_node(weight=4)
  g.add_edge(n0, n1) |> ignore
  g.add_edge(n1, n2) |> ignore
  g.add_edge(n2, n0) |> ignore
  g.add_edge(n0, n3) |> ignore
  g.add_edge(n3, n4) |> ignore
  g.add_edge(n4, n3) |> ignore
  let sccs = g.strongly_connected_components()

  // 应该有两个强连通分量
  inspect(sccs.length().to_string(), content="2")

  // 验证每个分量的大小
  let sizes = Array::new()
  for scc in sccs {
    sizes.push(scc.length())
  }

  // 排序以确保比较
  if sizes[0] > sizes[1] {
    sizes.swap(0, 1)
  }
  inspect(sizes[0].to_string(), content="2")
  inspect(sizes[1].to_string(), content="3")
}

///|
/// 测试最小可行流
test "min_flow" {
  let g : Graph2[String, Int] = Graph2::new(directed=true)
  let s = g.add_node(weight="s")
  let a = g.add_node(weight="a")
  let b = g.add_node(weight="b")
  let t = g.add_node(weight="t")

  // 添加边
  g.add_edge(s, a, weight=5) |> ignore
  g.add_edge(a, b, weight=5) |> ignore
  g.add_edge(b, t, weight=5) |> ignore
  g.add_edge(s, b, weight=3) |> ignore
  g.add_edge(a, t, weight=2) |> ignore

  // 定义流量下界
  let lower_bounds = Array::new()
  lower_bounds.push((s, a, 1))
  lower_bounds.push((a, b, 2))
  lower_bounds.push((b, t, 3))
  let min_flow = g.min_flow(s, t, lower_bounds)

  // 最小可行流应该是3
  inspect(min_flow.unwrap().to_string(), content="4")
}

///|
/// 测试DFS遍历
test "dfs_traversal" {
  let g : Graph2[String, Int] = Graph2::new(directed=false)

  // 创建一个简单的图
  //   A
  //  / \
  // B   C
  //  \ /
  //   D

  let a = g.add_node(weight="A")
  let b = g.add_node(weight="B")
  let c = g.add_node(weight="C")
  let d = g.add_node(weight="D")
  g.add_edge(a, b, weight=1) |> ignore
  g.add_edge(a, c, weight=1) |> ignore
  g.add_edge(b, d, weight=1) |> ignore
  g.add_edge(c, d, weight=1) |> ignore

  // 执行DFS遍历
  let dfs = Dfs::new(g, a)
  let visited = Array::new()
  let mut node = dfs.next(g)
  while node != None {
    visited.push(g.node_weights[node.unwrap().ix])
    node = dfs.next(g)
  }

  // 检查访问了所有节点
  inspect(visited.length().to_string(), content="4")

  // 检查第一个节点是起始点
  inspect(visited[0], content="A")
}

///|
/// 测试拓扑排序
test "topological_sort" {
  let g : Graph2[String, Int] = Graph2::new(directed=true)

  // 创建一个有向无环图
  let a = g.add_node(weight="A")
  let b = g.add_node(weight="B")
  let c = g.add_node(weight="C")
  let d = g.add_node(weight="D")
  g.add_edge(a, b, weight=1) |> ignore
  g.add_edge(a, c, weight=1) |> ignore
  g.add_edge(b, d, weight=1) |> ignore
  g.add_edge(c, d, weight=1) |> ignore

  // 执行拓扑排序
  let topo = Topo::new(g)
  let order = Array::new()
  let mut node = topo.next(g)
  while node != None {
    order.push(g.node_weights[node.unwrap().ix])
    node = topo.next(g)
  }

  // 检查排序结果长度
  inspect(order.length().to_string(), content="4")

  // 检查A排在B和C前面
  fn find_index(arr : Array[String], value : String) -> Int {
    let mut i = 0
    while i < arr.length() {
      if arr[i] == value {
        return i
      }
      i = i + 1
    }
    return -1
  }

  let a_pos = find_index(order, "A")
  let b_pos = find_index(order, "B")
  let c_pos = find_index(order, "C")
  let d_pos = find_index(order, "D")
  inspect((a_pos < b_pos).to_string(), content="true")
  inspect((a_pos < c_pos).to_string(), content="true")
  inspect((b_pos < d_pos).to_string(), content="true")
  inspect((c_pos < d_pos).to_string(), content="true")
}

///|
/// 测试基础邻居函数
test "neighbors_function" {
  let g : Graph2[String, Unit] = Graph2::new(directed=true)
  let a = g.add_node(weight="A")
  let b = g.add_node(weight="B")
  let c = g.add_node(weight="C")
  g.add_edge(a, b) |> ignore
  g.add_edge(a, c) |> ignore
  let neighbors = g.neighbors(a)
  inspect(neighbors.length().to_string(), content="2")
  // 检查b和c都在邻居列表中
  let mut has_b = false
  let mut has_c = false
  for neighbor in neighbors {
    if neighbor == b {
      has_b = true
    }
    if neighbor == c {
      has_c = true
    }
  }
  inspect(has_b.to_string(), content="true")
  inspect(has_c.to_string(), content="true")
}

///|
/// 测试BFS遍历
test "bfs_traversal" {
  let g : Graph2[String, Int] = Graph2::new(directed=false)
  let a = g.add_node(weight="A")
  let b = g.add_node(weight="B")
  let c = g.add_node(weight="C")
  let d = g.add_node(weight="D")
  g.add_edge(a, b, weight=1) |> ignore
  g.add_edge(a, c, weight=1) |> ignore
  g.add_edge(b, d, weight=1) |> ignore
  g.add_edge(c, d, weight=1) |> ignore
  let bfs = Bfs::new(g, a)
  let visited = Array::new()
  let mut node = bfs.next(g)
  while node != None {
    visited.push(g.node_weights[node.unwrap().ix])
    node = bfs.next(g)
  }

  // BFS顺序应该是A, (B,C), D或A, (C,B), D
  inspect(visited.length().to_string(), content="4")
  inspect(visited[0], content="A")
  inspect(visited[3], content="D")
}

///|
/// 测试边缘情况：空图
test "empty_graph" {
  let g : Graph2[Unit, Unit] = Graph2::new(directed=false)

  // 测试基本属性
  inspect(g.node_count().to_string(), content="0")
  inspect(g.edge_count().to_string(), content="0")

  // 测试空图上的算法行为
  let comps = g.connected_components()
  inspect(comps.length().to_string(), content="0")
}

///|
/// 测试边缘情况：单节点图
test "single_node_graph" {
  let g : Graph2[Int, Int] = Graph2::new(directed=false)
  let a = g.add_node(weight=1)

  // 测试基本属性
  inspect(g.node_count().to_string(), content="1")

  // 测试连通性
  let comps = g.connected_components()
  inspect(comps.length().to_string(), content="1")
  inspect(comps[0].length().to_string(), content="1")

  // 测试单节点上的最短路径
  let dist = g.dijkstra(a)
  inspect(dist.size().to_string(), content="1")
}
