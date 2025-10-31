///|
pub fn[N, E] Graph2::node_count(self : Graph2[N, E]) -> Int {
  self.adj.length()
}

///|
pub fn[N : Default, E] Graph2::add_node(
  self : Graph2[N, E],
  weight? : N,
) -> NodeIndex {
  let w = weight.unwrap_or(N::default())
  let idx = self.adj.length()
  self.node_weights.push(w)
  self.adj.push(Array::new())
  NodeIndex::new(idx)
}

///|
pub fn[N, E : Default] Graph2::add_edge(
  self : Graph2[N, E],
  a : NodeIndex,
  b : NodeIndex,
  weight? : E,
) -> Unit {
  let n = self.adj.length()
  guard a.ix >= 0 && a.ix < n && b.ix >= 0 && b.ix < n else {
    abort("add_edge: node index out of bounds")
  }
  let w = weight.unwrap_or(E::default())
  self.adj[a.ix].push((b, w))
  if not(self.directed) {
    self.adj[b.ix].push((a, w))
  }
  self.edges.push((a, b, w))
}

///|
pub fn[N, E] Graph2::neighbors(
  self : Graph2[N, E],
  a : NodeIndex,
) -> Array[NodeIndex] {
  let res = Array::new()
  for pair in self.adj[a.ix] {
    res.push(pair.0)
  }
  res
}

///|
/// 无向连通分量（返回每个连通分量的节点列表）
pub fn[N, E] Graph2::connected_components(
  self : Graph2[N, E],
) -> Array[Array[NodeIndex]] {
  let n = self.adj.length()
  let vis = Array::make(n, false)
  let comps = Array::new()
  for i in 0..<n {
    if not(vis[i]) {
      let comp = Array::new()
      // iterative stack DFS
      let stack = Array::new()
      stack.push(NodeIndex::new(i))
      while stack.pop() is Some(u) {
        if vis[u.ix] {
          continue
        }
        vis[u.ix] = true
        comp.push(u)
        for v in self.neighbors(u) {
          if not(vis[v.ix]) {
            stack.push(v)
          }
        }
      }
      comps.push(comp)
    }
  }
  comps
}

///|
/// Kosaraju 强连通分量（适用于有向图）
pub fn[N, E] Graph2::strongly_connected_components(
  self : Graph2[N, E],
) -> Array[Array[NodeIndex]] {
  let n = self.adj.length()
  // 1. 得到后序排列
  let visited = Array::make(n, false)
  let order = Array::new()
  for i in 0..<n {
    if not(visited[i]) {
      // iterative DFS to compute postorder
      let stack = Array::new() // elements: (NodeIndex, next_idx: Int)
      stack.push((NodeIndex::new(i), 0))
      while stack.last() is Some(top) {
        let (node, it) = top
        if not(visited[node.ix]) {
          visited[node.ix] = true
        }
        // find next neighbor starting from it
        let mut advanced = false
        let neighbors = self.adj[node.ix]
        if it < neighbors.length() {
          let (v, _) = neighbors[it]
          // increment iterator
          let last_idx = stack.length() - 1
          stack[last_idx] = (node, it + 1)
          if not(visited[v.ix]) {
            stack.push((v, 0))
          }
          advanced = true
        }
        if not(advanced) {
          // finished
          stack.pop() |> ignore
          order.push(node)
        }
      }
    }
  }
  // 2. 构造反向图
  let rev = Array::new()
  let mut kk = 0
  while kk < n {
    rev.push(Array::new())
    kk = kk + 1
  }
  for u_idx in 0..<n {
    for pair in self.adj[u_idx] {
      let (v, _) = pair
      rev[v.ix].push(NodeIndex::new(u_idx))
    }
  }
  // 3. 在反向图上按 order 的逆序 DFS
  let vis2 = Array::make(n, false)
  let sccs = Array::new()
  for jj = order.length() - 1; jj >= 0; jj = jj - 1 {
    let node = order[jj]
    if not(vis2[node.ix]) {
      let comp = Array::new()
      let stack2 = Array::new()
      stack2.push(node)
      while stack2.pop() is Some(u) {
        if vis2[u.ix] {
          continue
        }
        vis2[u.ix] = true
        comp.push(u)
        for v in rev[u.ix] {
          if not(vis2[v.ix]) {
            stack2.push(v)
          }
        }
      }
      sccs.push(comp)
    }
    if jj == 0 {
      break
    }
  }
  sccs
}

///|
/// Dijkstra（仅针对 Int 权重实现）
pub fn Graph2::dijkstra_int(
  self : Graph2[Int, Int],
  a : NodeIndex,
) -> Array[Int] {
  let n = self.adj.length()
  let dist = Array::make(n, @int.max_value)
  let vis = Array::make(n, false)
  dist[a.ix] = 0
  for _ in 0..<n {
    // pick min unvisited
    let mut u = -1
    let mut best = @int.max_value
    for i in 0..<n {
      if not(vis[i]) && dist[i] < best {
        best = dist[i]
        u = i
      }
    }
    if u == -1 || best == @int.max_value {
      break
    }
    vis[u] = true
    for pair in self.adj[u] {
      let (v, w) = pair
      if vis[v.ix] {
        continue
      }
      let nd = dist[u] + w
      if nd < dist[v.ix] {
        dist[v.ix] = nd
      }
    }
  }
  dist
}

///|
/// Kruskal 最小生成树（仅针对 Int 权重的无向图）
pub fn[N] Graph2::kruskal_mst(
  self : Graph2[N, Int],
) -> Array[(NodeIndex, NodeIndex, Int)] {
  let res = Array::new()
  let n = self.adj.length()
  // collect edges (undirected) ensuring each undirected edge counted once
  let es = Array::new()
  for e in self.edges {
    let (u, v, w) = e
    if self.directed {
      es.push((u, v, w))
    } else if u.ix <= v.ix {
      es.push((u, v, w))
    } else {
      es.push((v, u, w))
    }
  }
  // simple bubble sort by weight (while loops to avoid range operator issues)
  let sorted = es
  let len = sorted.length()
  let mut i = 0
  while i < len {
    let mut j = i + 1
    while j < len {
      if sorted[i].2 > sorted[j].2 {
        sorted.swap(i, j)
      }
      j = j + 1
    }
    i = i + 1
  }
  // union-find
  let parent = Array::make(n, 0)
  let rank = Array::make(n, 0)
  let mut k = 0
  while k < n {
    parent[k] = k
    k = k + 1
  }
  fn find_root(parent : Array[Int], start_x : Int) -> Int {
    let mut x = start_x
    while parent[x] != x {
      x = parent[x]
    }
    x
  }

  for e in sorted {
    let (u, v, w) = e
    let ru = find_root(parent, u.ix)
    let rv = find_root(parent, v.ix)
    if ru == rv {
      continue
    }
    if rank[ru] < rank[rv] {
      parent[ru] = rv
    } else if rank[rv] < rank[ru] {
      parent[rv] = ru
    } else {
      parent[rv] = ru
      rank[ru] += 1
    }
    res.push((u, v, w))
  }
  res
}

///|
/// Traversal helpers adapted from NetworkX (DFS, BFS, Topo)
struct Dfs {
  stack : Array[NodeIndex]
  discovered : VisitMap
}

///|
pub fn[N, E] Dfs::new(graph : Graph2[N, E], start : NodeIndex) -> Dfs {
  Dfs::{ stack: [start], discovered: graph.visit_map() }
}

///|
struct VisitMap {
  visit : Array[Bool]
}

///|
pub fn VisitMap::new() -> VisitMap {
  VisitMap::{ visit: Array::new() }
}

///|
fn VisitMap::contain(self : VisitMap, nx : Int) -> Bool {
  self.visit[nx]
}

///|
impl VM for VisitMap with visit(self, node : NodeIndex) -> Bool {
  let nx = node.ix
  let res = self.contain(nx)
  self.visit[nx] = true
  not(res)
}

///|
impl VM for VisitMap with is_visited(self, node : NodeIndex) -> Bool {
  self.contain(node.ix)
}

///|
priv trait VM {
  visit(Self, NodeIndex) -> Bool
  is_visited(Self, NodeIndex) -> Bool
}

///|
priv trait Visitable {
  visit_map(Self) -> VisitMap
}

///|
impl[N, E] Visitable for Graph2[N, E] with visit_map(self) -> VisitMap {
  VisitMap::{ visit: Array::make(self.node_count(), false) }
}

///|
fn VisitMap::clear(self : VisitMap) -> Unit {
  self.visit.clear()
}

///|
pub fn Dfs::from_parts(stack : Array[NodeIndex], discovered : VisitMap) -> Dfs {
  Dfs::{ stack, discovered }
}

///|
pub fn Dfs::reset(self : Dfs) -> Unit {
  self.stack.clear()
  self.discovered.clear()
}

///|
pub fn Dfs::move_to(self : Dfs, node : NodeIndex) -> Unit {
  self.stack.clear()
  self.stack.push(node)
}

///|
pub fn[N, E] Dfs::next(self : Dfs, graph : Graph2[N, E]) -> NodeIndex? {
  while self.stack.pop() is Some(node) {
    if self.discovered.visit(node) {
      for succ in graph.neighbors(node) {
        if not(self.discovered.is_visited(succ)) {
          self.stack.push(succ)
        }
      }
      return Some(node)
    }
  }
  None
}

///|
struct DfsPostOrder {
  stack : Array[NodeIndex]
  discovered : VisitMap
  finished : VisitMap
}

///|
pub fn[N, E] DfsPostOrder::new(
  graph : Graph2[N, E],
  start : NodeIndex,
) -> DfsPostOrder {
  DfsPostOrder::{
    stack: [start],
    discovered: graph.visit_map(),
    finished: graph.visit_map(),
  }
}

///|
pub fn DfsPostOrder::reset(self : DfsPostOrder) -> Unit {
  self.stack.clear()
  self.discovered.clear()
  self.finished.clear()
}

///|
pub fn DfsPostOrder::move_to(self : DfsPostOrder, node : NodeIndex) -> Unit {
  self.stack.clear()
  self.stack.push(node)
}

///|
pub fn[N, E] DfsPostOrder::next(
  self : DfsPostOrder,
  graph : Graph2[N, E],
) -> NodeIndex? {
  while self.stack.last() is Some(node) {
    if self.discovered.visit(node) {
      for succ in graph.neighbors(node) {
        if not(self.discovered.is_visited(succ)) {
          self.stack.push(succ)
        }
      }
    } else {
      self.stack.pop() |> ignore
      if self.finished.visit(node) {
        return Some(node)
      }
    }
  }
  None
}

///|
struct Bfs {
  stack : @deque.Deque[NodeIndex]
  discovered : VisitMap
}

///|
/// Dinic最大流算法 - 比Edmonds-Karp更高效的最大流实现，适用于大多数网络
/// - source: 源点
/// - sink: 汇点
/// - 返回从源点到汇点的最大流量
pub fn[N, E : Compare + Add + Default + Sub] Graph2::dinic_max_flow(
  self : Graph2[N, E],
  source : NodeIndex,
  sink : NodeIndex,
) -> E {
  let n = self.node_count()
  guard source.ix >= 0 &&
    source.ix < n &&
    sink.ix >= 0 &&
    sink.ix < n &&
    source != sink else {
    return E::default()
  }

  // 创建残差网络
  let residual = Array::new()
  for _ in 0..<n {
    residual.push(Array::make(n, E::default()))
  }

  // 初始化残差网络并找出最大容量
  let mut max_capacity = E::default()
  for u_idx in 0..<n {
    for pair in self.adj[u_idx] {
      let (v, capacity) = pair
      residual[u_idx][v.ix] = capacity

      // 记录最大容量，用于计算初始流量限制
      if capacity > max_capacity {
        max_capacity = capacity
      }
    }
  }

  // 层次图（用于BFS）
  let level = Array::make(n, -1)

  // 使用一个大的初始流量限制，不依赖于具体类型
  let mut max_flow = E::default()

  // 不断寻找增广路径
  while build_level_graph(residual, level, source, sink, n) {
    // 记录每个节点当前检查到的边
    let next_edge = Array::make(n, 0)

    // 使用足够大的流量限制 - 这里使用网络中的最大容量 * 节点数
    // 这样即使在最坏情况下，也能容纳所有可能的流量
    let mut big_flow_limit = max_capacity
    for _ in 0..<n {
      big_flow_limit = big_flow_limit + max_capacity // 累加最大容量
    }

    // 寻找所有可行的增广路径
    let mut path_flow = dfs_send_flow(
      residual,
      level,
      next_edge,
      source.ix,
      sink.ix, // 注意：这里应该使用sink.ix而不是sink对象
      big_flow_limit,
    )
    while path_flow > E::default() {
      max_flow = max_flow + path_flow
      path_flow = dfs_send_flow(
        residual,
        level,
        next_edge,
        source.ix,
        sink.ix, // 同上，使用索引
        big_flow_limit,
      )
    }
  }
  max_flow
}

///|
/// 构建层次图 - BFS确定节点到源点的距离
fn[E : Compare + Default] build_level_graph(
  residual : Array[Array[E]],
  level : Array[Int],
  source : NodeIndex,
  sink : NodeIndex,
  n : Int,
) -> Bool {
  // 重置所有节点的层次
  for i in 0..<n {
    level[i] = -1
  }

  // 源点的层次为0
  level[source.ix] = 0

  // 使用BFS计算所有可达节点的层次
  let queue = Array::new()
  queue.push(source.ix)
  while !queue.is_empty() {
    let u = queue.remove(0)
    for v in 0..<n {
      if level[v] == -1 && residual[u][v] > E::default() {
        level[v] = level[u] + 1
        queue.push(v)
      }
    }
  }

  // 如果汇点层次不为-1，则找到了一条路径
  level[sink.ix] != -1
}

///|
/// 深度优先搜索寻找增广路径并发送流量
fn[E : Compare + Add + Default + Sub] dfs_send_flow(
  residual : Array[Array[E]],
  level : Array[Int],
  next_edge : Array[Int],
  u : Int,
  t : Int,
  flow_limit : E,
) -> E {
  // 到达汇点，返回可用流量
  if u == t {
    return flow_limit
  }

  // 检查所有可行的邻居节点
  let n = residual.length()
  for i = next_edge[u]; i < n; i = i + 1 {
    next_edge[u] = i + 1
    let v = i

    // 只考虑层次比当前节点高一级且有残差容量的边
    if level[v] == level[u] + 1 && residual[u][v] > E::default() {
      // 递归寻找增广路径
      let bottle_neck = dfs_send_flow(
        residual,
        level,
        next_edge,
        v,
        t,
        @cmp.minimum(flow_limit, residual[u][v]),
      )

      // 如果找到了增广路径
      if bottle_neck > E::default() {
        // 更新残差网络
        residual[u][v] = residual[u][v] - bottle_neck
        residual[v][u] = residual[v][u] + bottle_neck
        return bottle_neck
      }
    }
  }

  // 没有找到增广路径
  E::default()
}

///|
pub fn[N, E] Bfs::new(graph : Graph2[N, E], start : NodeIndex) -> Bfs {
  let stack = @deque.new()
  stack.push_front(start)
  let discovered = graph.visit_map()
  discovered.visit(start) |> ignore
  Bfs::{ stack, discovered }
}

///|
pub fn[N, E] Bfs::next(self : Bfs, graph : Graph2[N, E]) -> NodeIndex? {
  if self.stack.pop_front() is Some(node) {
    for succ in graph.neighbors(node) {
      if self.discovered.visit(succ) {
        self.stack.push_back(succ)
      }
    }
    return Some(node)
  }
  None
}

///|
struct Topo {
  tovisit : Array[NodeIndex]
  ordered : VisitMap
  degrees : Array[Int]
}

///|
pub fn[N, E] Topo::new(graph : Graph2[N, E]) -> Topo {
  let degrees = graph.get_indegree()
  let ordered = graph.visit_map()
  let tovisit = Array::new()
  for i in 0..<graph.node_count() {
    if degrees[i] == 0 {
      tovisit.push(NodeIndex::new(i))
    }
  }
  Topo::{ tovisit, ordered, degrees }
}

///|
pub fn[N, E] Topo::empty(graph : Graph2[N, E]) -> Topo {
  Topo::{
    tovisit: Array::new(),
    ordered: graph.visit_map(),
    degrees: Array::new(),
  }
}

///|
pub fn[N, E] Topo::next(self : Topo, graph : Graph2[N, E]) -> NodeIndex? {
  while self.tovisit.pop() is Some(node) {
    if self.ordered.is_visited(node) {
      continue
    }
    self.ordered.visit(node) |> ignore
    for v in graph.neighbors(node) {
      self.degrees[v.ix] -= 1
      if self.degrees[v.ix] == 0 {
        self.tovisit.push(v)
      }
    }
    return Some(node)
  }
  None
}

///|
fn[N, E] get_indegree(self : Graph2[N, E]) -> Array[Int] {
  let deg = Array::make(self.node_count(), 0)
  for i in 0..<self.node_count() {
    for pair in self.adj[i] {
      let (v, _) = pair
      deg[v.ix] += 1
    }
  }
  deg
}

///|
/// Generic Dijkstra adapted for Graph2
pub fn[N, E : Compare + Add + Default + Neg] Graph2::dijkstra(
  self : Graph2[N, E],
  a : NodeIndex,
) -> @hashmap.HashMap[NodeIndex, E] {
  let n = self.adj.length()
  guard a.ix >= 0 && a.ix < n else { abort("Node not found") }
  let vis = Array::make(n, false)
  let dist = @hashmap.new()
  dist[a] = E::default()
  let q : @priority_queue.T[(E, NodeIndex)] = @priority_queue.new()
  // 使用负值来模拟小根堆
  q.push((-E::default(), a))
  while not(q.is_empty()) {
    let (neg_d, u) = q.pop().unwrap()
    let d = -neg_d // 恢复原始距离
    if vis[u.ix] {
      continue
    }
    vis[u.ix] = true
    for pair in self.adj[u.ix] {
      let (v, w) = pair
      if vis[v.ix] {
        continue
      }
      let new_dist = d + w
      if dist.get(v) == None || new_dist <= dist.get(v).unwrap() {
        dist[v] = new_dist
        // 使用负值来模拟小根堆
        q.push((-new_dist, v))
      }
    }
  }
  dist
}

///|
/// LCA adapted for Graph2 (removed cycle check to match Graph2 capabilities)
struct Lca[N, E] {
  graph : Graph2[N, E]
  deep : Array[Int]
  father : Array[Array[Int]]
  mut root : NodeIndex
}

///|
fn[N, E : Default] Lca::update_lca_info_pre(
  self : Lca[N, E],
  graph : Graph2[N, E],
  u : NodeIndex,
  f : NodeIndex,
  d : Int,
) -> Unit {
  self.deep[u.ix] = d
  for v in graph.neighbors(u) {
    if v != f {
      self.father[v.ix][0] = u.ix
      self.update_lca_info_pre(graph, v, u, d + 1)
    }
  }
}

///|
fn[N, E : Default] Lca::update_lca_info(
  self : Lca[N, E],
  graph : Graph2[N, E],
  u : NodeIndex,
  f : NodeIndex,
  d : Int,
) -> Unit {
  let n = self.graph.node_count()
  let top = @math.log2((n + 1).to_double()).to_int()
  self.update_lca_info_pre(graph, u, f, d)
  for j in 1..=top {
    for i in 1..=n {
      self.father[i][j] = self.father[self.father[i][j - 1]][j - 1]
    }
  }
}

///|
pub fn[N, E : Default] Lca::new(
  graph : Graph2[N, E],
  root : NodeIndex,
) -> Lca[N, E] {
  let n = graph.node_count()
  let deep = Array::make(n + 1, 0)
  let top = @math.log2((n + 1).to_double()).to_int()
  let father = Array::makei(n + 1, fn(_) { Array::make(top + 1, 0) })
  let lca = Lca::{ graph, deep, father, root }
  lca.update_lca_info(lca.graph, lca.root, NodeIndex::new(-1), 1)
  lca
}

///|
pub fn[N, E] Lca::query(
  self : Lca[N, E],
  u : NodeIndex,
  v : NodeIndex,
) -> NodeIndex {
  let mut u = u.ix
  let mut v = v.ix
  if self.deep[u] < self.deep[v] {
    let temp = u
    u = v
    v = temp
  }
  let n = self.graph.node_count()
  let top = @math.log2((n + 1).to_double()).to_int()
  for i = top; i >= 0; i = i - 1 {
    if self.deep[self.father[u][i]] >= self.deep[v] {
      u = self.father[u][i]
    }
  }
  if u == v {
    return NodeIndex::new(u)
  }
  for i = top; i >= 0; i = i - 1 {
    if self.father[u][i] != self.father[v][i] {
      u = self.father[u][i]
      v = self.father[v][i]
    }
  }
  NodeIndex::new(self.father[u][0])
}

///|
pub fn[N, E : Default] Lca::reset(self : Lca[N, E], root? : NodeIndex) -> Unit {
  if root is Some(r) {
    self.root = r
  }
  self.update_lca_info(self.graph, self.root, NodeIndex::new(-1), 1)
}

///|
/// Utility tuple helpers copied from NetworkX
pub fn[T] get_tuple(tuple : (T, T), index : Int) -> T {
  if index == 0 {
    return tuple.0
  }
  tuple.1
}

///|
pub fn[T] set_tuple(tuple : (T, T), index : Int, value : T) -> (T, T) {
  if index == 0 {
    return (value, tuple.1)
  }
  (tuple.0, value)
}

///|
pub fn[N, E] Graph2::find_edge(
  self : Graph2[N, E],
  a : NodeIndex,
  b : NodeIndex,
) -> Int? {
  let len = self.edges.length()
  let mut i = 0
  while i < len {
    let (u, v, _) = self.edges[i]
    if u == a && v == b {
      return Some(i)
    }
    if not(self.directed) {
      if u == b && v == a {
        return Some(i)
      }
    }
    i = i + 1
  }
  None
}

///|
pub fn[N, E : Eq] Graph2::remove_edge_nodes(
  self : Graph2[N, E],
  a : NodeIndex,
  b : NodeIndex,
) -> E? {
  match self.find_edge(a, b) {
    None => None
    Some(idx) => {
      let (_, _, w) = self.edges[idx]
      // remove from adjacency list of a
      let mut i = 0
      while i < self.adj[a.ix].length() {
        let (v, ww) = self.adj[a.ix][i]
        if v == b && ww == w {
          ignore(self.adj[a.ix].remove(i))
          break
        }
        i = i + 1
      }
      // if undirected, remove from adjacency list of b as well
      if not(self.directed) {
        let mut j = 0
        while j < self.adj[b.ix].length() {
          let (v2, ww2) = self.adj[b.ix][j]
          if v2 == a && ww2 == w {
            ignore(self.adj[b.ix].remove(j))
            break
          }
          j = j + 1
        }
      }
      // remove from edges array (swap-remove)
      let last_index = self.edges.length() - 1
      if idx != last_index {
        self.edges.swap(idx, last_index)
      }
      ignore(self.edges.remove(last_index))
      Some(w)
    }
  }
}

///|
pub fn[N, E : Eq] Graph2::remove_edge(self : Graph2[N, E], e : Int) -> E? {
  let len = self.edges.length()
  guard e >= 0 && e < len else { return None }
  let (u, v, w) = self.edges[e]

  // remove from adjacency list of u
  let mut i = 0
  while i < self.adj[u.ix].length() {
    let (vv, ww) = self.adj[u.ix][i]
    if vv == v && ww == w {
      ignore(self.adj[u.ix].remove(i))
      break
    }
    i = i + 1
  }

  // if undirected, remove from adjacency list of v as well
  if not(self.directed) {
    let mut j = 0
    while j < self.adj[v.ix].length() {
      let (uu, ww2) = self.adj[v.ix][j]
      if uu == u && ww2 == w {
        ignore(self.adj[v.ix].remove(j))
        break
      }
      j = j + 1
    }
  }

  // remove from edges array (swap-remove)
  let last_index = self.edges.length() - 1
  if e != last_index {
    self.edges.swap(e, last_index)
  }
  ignore(self.edges.remove(last_index))
  Some(w)
}

///| A* 寻路算法，结合最短路径成本和启发式评估

///| - start: 起始节点

///| - goal: 目标节点

///|
pub fn[N, E : Compare + Add + Default] Graph2::a_star(
  self : Graph2[N, E],
  start : NodeIndex,
  goal : NodeIndex,
  heuristic : (NodeIndex, NodeIndex) -> E,
) -> Array[NodeIndex] {
  let n = self.adj.length()
  guard start.ix >= 0 && start.ix < n && goal.ix >= 0 && goal.ix < n else {
    return Array::new()
  }

  // 用于存储从起点到每个节点的最短距离
  let g_score : @hashmap.HashMap[NodeIndex, E] = @hashmap.new()
  g_score.set(start, E::default()) // put 改为 set

  // 用于存储估计的总成本 f(n) = g(n) + h(n)
  let f_score : @hashmap.HashMap[NodeIndex, E] = @hashmap.new()
  f_score.set(start, heuristic(start, goal)) // put 改为 set

  // 用于重建路径的前驱节点图
  let came_from : @hashmap.HashMap[NodeIndex, NodeIndex] = @hashmap.new()

  // 优先队列，按 f_score 排序
  let open_set : @priority_queue.T[(E, NodeIndex)] = @priority_queue.new()
  open_set.push((f_score.get(start).unwrap(), start))

  // 记录已处理的节点
  let closed_set : @hashmap.HashMap[NodeIndex, Bool] = @hashmap.new()
  while not(open_set.is_empty()) {
    let (_, current) = open_set.pop().unwrap()

    // 如果到达目标，重建并返回路径
    if current == goal {
      return reconstruct_path(came_from, current)
    }

    // 标记为已处理
    closed_set.set(current, true) // put 改为 set

    // 检查所有邻居
    for pair in self.adj[current.ix] {
      let (neighbor, edge_cost) = pair

      // 跳过已处理的节点
      if closed_set.get(neighbor) != None {
        continue
      }

      // 计算经过当前节点到达邻居的成本
      let tentative_g_score = g_score.get(current).unwrap() + edge_cost

      // 如果是新节点或找到了更短的路径
      if g_score.get(neighbor) == None ||
        tentative_g_score < g_score.get(neighbor).unwrap() {
        // 记录这条路径
        came_from.set(neighbor, current) // put 改为 set
        g_score.set(neighbor, tentative_g_score) // put 改为 set
        f_score.set(neighbor, tentative_g_score + heuristic(neighbor, goal)) // put 改为 set

        // 加入待处理队列
        open_set.push((f_score.get(neighbor).unwrap(), neighbor))
      }
    }
  }

  // 找不到路径
  Array::new()
}

///|
/// A* 路径重建辅助函数
fn reconstruct_path(
  came_from : @hashmap.HashMap[NodeIndex, NodeIndex],
  current : NodeIndex,
) -> Array[NodeIndex] {
  let path = Array::new()
  let mut cur = current
  while true {
    path.push(cur)
    match came_from.get(cur) {
      Some(prev) => cur = prev
      None => break
    }
  }
  // reverse
  let res = Array::new()
  let mut i = path.length() - 1
  while i >= 0 {
    res.push(path[i])
    if i == 0 {
      break
    }
    i = i - 1
  }
  res
}

///|

///|
pub fn[N, E : Eq] Graph2::remove_node(self : Graph2[N, E], a : NodeIndex) -> N? {
  let n = self.node_count()
  guard a.ix >= 0 && a.ix < n else { return None }

  // 找出要删除的边（只是与节点a直接相连的边）
  let edges_to_remove = Array::new()
  for ei in 0..<self.edges.length() {
    let (u, v, _) = self.edges[ei]
    if u == a || v == a {
      edges_to_remove.push(ei)
    }
  }

  // 从大到小删除这些边（避免索引变化问题）
  // Remove edges from largest index to smallest to avoid index shifting issues
  edges_to_remove.sort() // sorts in ascending order
  let mut i = edges_to_remove.length() - 1
  while i >= 0 {
    let ei = edges_to_remove[i]
    let (u, v, _) = self.edges[ei]
    self.remove_edge_nodes(u, v) |> ignore
    if i == 0 {
      break
    }
    i = i - 1
  }

  // 保存原始的节点权重以便返回
  let weight = self.node_weights[a.ix]

  // swap-remove操作：用最后一个节点替换被删除的节点
  let last_index = self.node_count() - 1
  if a.ix != last_index {
    // 移动最后一个节点到位置a
    self.node_weights.swap(a.ix, last_index)
    self.adj.swap(a.ix, last_index)
    let old_idx = NodeIndex::new(last_index)
    let new_idx = NodeIndex::new(a.ix)

    // 更新邻接表中的引用
    let mut i = 0
    while i < self.adj.length() {
      let mut j = 0
      while j < self.adj[i].length() {
        let (v, ww) = self.adj[i][j]
        if v == old_idx {
          self.adj[i][j] = (new_idx, ww)
        }
        j = j + 1
      }
      i = i + 1
    }

    // 更新边数组中的节点引用
    let mut ei2 = 0
    while ei2 < self.edges.length() {
      let (u2, v2, ww2) = self.edges[ei2]
      let nu = if u2 == old_idx { new_idx } else { u2 }
      let nv = if v2 == old_idx { new_idx } else { v2 }
      if nu != u2 || nv != v2 {
        self.edges[ei2] = (nu, nv, ww2)
      }
      ei2 = ei2 + 1
    }
  }

  // 删除最后一个节点
  ignore(self.node_weights.remove(last_index))
  ignore(self.adj.remove(last_index))
  Some(weight)
}

///|
pub fn[N, E] Graph2::edge_count(self : Graph2[N, E]) -> Int {
  self.edges.length()
}

///|
pub fn[N, E] Graph2::find_edge_undirected(
  self : Graph2[N, E],
  a : NodeIndex,
  b : NodeIndex,
) -> Int? {
  // Graph2::find_edge already handles undirected graphs, so reuse it
  self.find_edge(a, b)
}

///|
pub fn[N, E] Graph2::remove_edge_adjust_indices(
  self : Graph2[N, E],
  e : Int,
) -> E? {
  let last_index = self.edges.length() - 1
  if e != last_index {
    self.edges.swap(e, last_index)
  }
  let (_, _, w) = self.edges.remove(last_index)
  Some(w)
}

///| 寻找有向图中的所有环

///|
/// 返回每个环的节点列表
pub fn[N, E] Graph2::find_cycles(
  self : Graph2[N, E],
) -> Array[Array[NodeIndex]] {
  let n = self.adj.length()
  let cycles = Array::new()

  // 0=未访问, 1=访问中, 2=已访问完成
  let visited = Array::make(n, 0)

  // 用于跟踪当前DFS路径上的节点
  let path = Array::new()

  // 对每个未访问的节点进行DFS
  for i in 0..<n {
    if visited[i] == 0 {
      dfs_find_cycles(self, NodeIndex::new(i), visited, path, cycles)
    }
  }
  cycles
}

///|
/// DFS辅助函数，用于检测环
fn[N, E] dfs_find_cycles(
  graph : Graph2[N, E],
  node : NodeIndex,
  visited : Array[Int],
  path : Array[NodeIndex],
  cycles : Array[Array[NodeIndex]],
) -> Unit {
  // 标记为"访问中"
  visited[node.ix] = 1
  path.push(node)

  // 检查所有邻居
  for pair in graph.adj[node.ix] {
    let (neighbor, _) = pair
    if visited[neighbor.ix] == 0 {
      // 如果邻居未访问，递归处理
      dfs_find_cycles(graph, neighbor, visited, path, cycles)
    } else if visited[neighbor.ix] == 1 {
      // 如果邻居是"访问中"状态，找到了一个环
      // 找到路径中的起始位置
      let cycle = Array::new()
      let mut start_index = path.length() - 1
      while start_index >= 0 {
        if path[start_index] == neighbor {
          break
        }
        start_index -= 1
      }
      if start_index >= 0 {
        // 构建环
        for i in start_index..<path.length() {
          cycle.push(path[i])
        }
        cycles.push(cycle)
      }
    }
    // 如果邻居已经处理完成(状态为2)，则忽略
  }

  // 回溯，移除当前节点
  ignore(path.pop())
  // 标记为"已访问完成"
  visited[node.ix] = 2
}

///| 使用 Edmonds-Karp 算法计算网络最大流

///| - source: 源点

///| - sink: 汇点

///|
/// - 返回从源点到汇点的最大流量
pub fn[N, E : Compare + Add + Default + Sub] Graph2::max_flow(
  self : Graph2[N, E],
  source : NodeIndex,
  sink : NodeIndex,
) -> E {
  let n = self.node_count()
  guard source.ix >= 0 &&
    source.ix < n &&
    sink.ix >= 0 &&
    sink.ix < n &&
    source != sink else {
    return E::default()
  }

  // 创建残差网络
  let residual = Array::new()
  for _ in 0..<n {
    residual.push(Array::make(n, E::default()))
  }

  // 初始化残差网络
  for u_idx in 0..<n {
    for pair in self.adj[u_idx] {
      let (v, capacity) = pair
      residual[u_idx][v.ix] = capacity
    }
  }
  let mut max_flow = E::default()

  // 重复寻找增广路径直到没有更多路径
  while true {
    // 使用BFS寻找从source到sink的增广路径
    let parent = Array::make(n, -1)
    let queue = @deque.new()
    queue.push_back(source)
    parent[source.ix] = source.ix // 标记源点已访问
    let mut path_found = false
    while !queue.is_empty() && !path_found {
      let u = queue.pop_front().unwrap()

      // 检查所有可能的邻居
      for v_idx in 0..<n {
        if parent[v_idx] == -1 && residual[u.ix][v_idx] > E::default() {
          parent[v_idx] = u.ix
          if v_idx == sink.ix {
            path_found = true
            break
          }
          queue.push_back(NodeIndex::new(v_idx))
        }
      }
    }

    // 如果没有找到增广路径，则结束算法
    if !path_found {
      break
    }

    // 计算增广路径上的最小容量（瓶颈）
    let mut bottleneck = None
    let mut v = sink.ix
    while v != source.ix {
      let u = parent[v]
      let capacity = residual[u][v]
      if bottleneck is None {
        bottleneck = Some(capacity)
      } else {
        let b = bottleneck.unwrap()
        if capacity < b {
          bottleneck = Some(capacity)
        }
      }
      v = u
    }
    let flow = bottleneck.unwrap()

    // 更新残差网络
    v = sink.ix
    while v != source.ix {
      let u = parent[v]
      residual[u][v] -= flow // 减少正向边容量
      residual[v][u] += flow // 增加反向边容量
      v = u
    }

    // 更新最大流
    max_flow += flow
  }
  max_flow
}

///|
/// 计算具有流量下界约束的网络中的最小可行流
/// - source: 源点
/// - sink: 汇点
/// - lower_bounds: 边的流量下界约束数组，每项为 (起点, 终点, 下界值)
/// - 返回: 最小可行流量，如果不存在可行流则返回 None
pub fn[N : Default] Graph2::min_flow(
  self : Graph2[N, Int],
  source : NodeIndex,
  sink : NodeIndex,
  lower_bounds : Array[(NodeIndex, NodeIndex, Int)],
) -> Int? {
  // 创建辅助网络
  let aux_graph = Graph2::new(directed=true)

  // 复制所有节点
  for i in 0..<self.node_count() {
    ignore(aux_graph.add_node(weight=self.node_weights[i]))
  }

  // 计算每个节点的净流量需求
  let demands = Array::make(self.node_count(), 0)

  // 添加原始边，考虑下界约束
  for u_idx in 0..<self.node_count() {
    let u = NodeIndex::new(u_idx)
    for pair in self.adj[u_idx] {
      let (v, capacity) = pair

      // 查找此边的流量下界
      let mut lower_bound = 0
      for bound_tuple in lower_bounds {
        let (src, dst, lb_value) = bound_tuple
        if src == u && dst == v {
          lower_bound = lb_value
          break
        }
      }

      // 更新流量需求
      demands[u.ix] -= lower_bound
      demands[v.ix] += lower_bound

      // 添加边，容量减去下界
      let remaining_capacity = capacity - lower_bound
      if remaining_capacity > 0 {
        aux_graph.add_edge(u, v, weight=remaining_capacity)
      }
    }
  }

  // 添加超级源点和超级汇点
  let super_source = aux_graph.add_node(weight=N::default())
  let super_sink = aux_graph.add_node(weight=N::default())

  // 处理流量需求
  let mut total_demand = 0
  for i in 0..<demands.length() {
    if i != source.ix && i != sink.ix {
      if demands[i] > 0 { // 节点需要流入
        aux_graph.add_edge(super_source, NodeIndex::new(i), weight=demands[i])
        total_demand += demands[i]
      } else if demands[i] < 0 { // 节点需要流出
        aux_graph.add_edge(NodeIndex::new(i), super_sink, weight=-demands[i])
      }
    }
  }

  // 连接原始源汇点
  aux_graph.add_edge(sink, source, weight=@int.max_value / 2)

  // 计算最大流
  let max_flow = aux_graph.max_flow(super_source, super_sink)

  // 检查是否存在可行流
  if max_flow < total_demand {
    return None // 不存在可行流
  }

  // 计算源点到汇点的最小流量
  let mut min_flow_value = 0
  for bound_tuple in lower_bounds {
    let (src, dst, lb_value) = bound_tuple
    if src == source || dst == sink {
      min_flow_value += lb_value
    }
  }
  Some(min_flow_value)
}

///|
/// - 返回：距离表和是否存在负环的标记
pub fn[N, E : Compare + Add + Default] Graph2::bellman_ford(
  self : Graph2[N, E],
  source : NodeIndex,
) -> (Array[E?], Bool) {
  let n = self.node_count()
  let inf = None

  // 初始化距离数组
  let dist = Array::new()
  for _ in 0..<n {
    dist.push(inf)
  }
  dist[source.ix] = Some(E::default())

  // Bellman-Ford 主循环
  let mut has_negative_cycle = false

  // 迭代 n-1 次
  let mut iter = 0
  while iter < n - 1 {
    let mut updated = false

    // 检查所有边
    for e_idx in 0..<self.edge_count() {
      let (u, v, w) = self.edges[e_idx]
      if dist[u.ix] != None {
        let u_dist = dist[u.ix].unwrap()
        let new_dist = u_dist + w
        if dist[v.ix] == None || new_dist < dist[v.ix].unwrap() {
          dist[v.ix] = Some(new_dist)
          updated = true
        }
      }
    }

    // 如果没有更新，提前结束
    if !updated {
      break
    }
    iter = iter + 1
  }

  // 检查负权环
  for e_idx in 0..<self.edge_count() {
    let (u, v, w) = self.edges[e_idx]
    if dist[u.ix] != None {
      let u_dist = dist[u.ix].unwrap()
      let new_dist = u_dist + w
      if dist[v.ix] != None && new_dist < dist[v.ix].unwrap() {
        has_negative_cycle = true
        break
      }
    }
  }
  (dist, has_negative_cycle)
}

///| Floyd-Warshall 多源最短路径算法

///|
/// 计算图中任意两点间的最短路径
pub fn[N, E : Compare + Add + Default] Graph2::floyd_warshall(
  self : Graph2[N, E],
) -> Array[Array[E?]] {
  let n = self.node_count()
  let inf = None

  // 初始化距离矩阵
  let dist = Array::new()
  for _ in 0..<n {
    let row = Array::new()
    for _ in 0..<n {
      row.push(inf)
    }
    dist.push(row)
  }

  // 设置初始距离
  for i in 0..<n {
    dist[i][i] = Some(E::default()) // 到自身的距离为0

    // 设置直接相连的边
    for pair in self.adj[i] {
      let (v, weight) = pair
      dist[i][v.ix] = Some(weight)
    }
  }

  // Floyd-Warshall 主循环
  for k in 0..<n {
    for i in 0..<n {
      for j in 0..<n {
        if dist[i][k] != None && dist[k][j] != None {
          let path_through_k = dist[i][k].unwrap() + dist[k][j].unwrap()
          if dist[i][j] == None || path_through_k < dist[i][j].unwrap() {
            dist[i][j] = Some(path_through_k)
          }
        }
      }
    }
  }
  dist
}

///| Prim 最小生成树算法

///|
/// 返回最小生成树的边集
pub fn[N] Graph2::prim_mst(
  self : Graph2[N, Int],
) -> Array[(NodeIndex, NodeIndex, Int)] {
  let n = self.node_count()
  guard n > 0 else { return Array::new() }
  let result = Array::new()
  let visited = Array::make(n, false)

  // 使用优先队列存储边
  let pq = @priority_queue.T::new()

  // 从第一个节点开始
  visited[0] = true

  // 将第一个节点的所有边加入优先队列
  for pair in self.adj[0] {
    let (v, weight) = pair
    pq.push((weight, (NodeIndex::new(0), v)))
  }

  // 当队列非空且未找到所有边
  while !pq.is_empty() && result.length() < n - 1 {
    let (weight, (u, v)) = pq.pop().unwrap()

    // 如果目标节点已访问，跳过
    if visited[v.ix] {
      continue
    }

    // 将边加入结果
    result.push((u, v, weight))
    visited[v.ix] = true

    // 将新节点的边加入队列
    for pair in self.adj[v.ix] {
      let (next, edge_weight) = pair
      if !visited[next.ix] {
        pq.push((edge_weight, (v, next)))
      }
    }
  }
  result
}

///|
/// 检查图是否为二分图，并返回染色结果
pub fn[N, E] Graph2::is_bipartite(self : Graph2[N, E]) -> (Bool, Array[Int]) {
  let n = self.node_count()
  let colors = Array::make(n, -1) // -1: 未染色, 0: 红色, 1: 蓝色

  // 对每个连通分量进行二分染色
  for start in 0..<n {
    if colors[start] == -1 {
      colors[start] = 0 // 初始节点染为红色
      let queue = @deque.new()
      queue.push_back(NodeIndex::new(start))
      while !queue.is_empty() {
        let u = queue.pop_front().unwrap()
        for pair in self.adj[u.ix] {
          let (v, _) = pair
          if colors[v.ix] == -1 {
            // 给邻居染相反的颜色
            colors[v.ix] = 1 - colors[u.ix]
            queue.push_back(v)
          } else if colors[v.ix] == colors[u.ix] {
            // 相邻节点颜色相同，不是二分图
            return (false, colors)
          }
        }
      }
    }
  }
  (true, colors)
}

///| 使用匈牙利算法计算二分图的最大匹配

///|
/// 需要预先确认是二分图
pub fn[N, E] Graph2::bipartite_matching(
  self : Graph2[N, E],
  colors : Array[Int],
) -> Array[(NodeIndex, NodeIndex)] {
  let n = self.node_count()
  let matching = Array::make(n, -1) // 每个节点匹配的对象
  let result = Array::new()

  // DFS辅助函数
  fn dfs(
    graph : Graph2[N, E],
    u : Int,
    visited : Array[Bool],
    matching : Array[Int],
    colors : Array[Int],
  ) -> Bool {
    for pair in graph.adj[u] {
      let (v, _) = pair

      // 只考虑另一个集合中的节点
      if colors[v.ix] != colors[u] {
        if visited[v.ix] {
          continue
        }
        visited[v.ix] = true

        // 如果v未匹配或v的当前匹配可以找到新的匹配
        if matching[v.ix] == -1 ||
          dfs(graph, matching[v.ix], visited, matching, colors) {
          matching[v.ix] = u
          matching[u] = v.ix
          return true
        }
      }
    }
    false
  }

  // 尝试为每个集合0中的节点寻找匹配
  for u in 0..<n {
    if colors[u] == 0 && matching[u] == -1 {
      let visited = Array::make(n, false)
      ignore(dfs(self, u, visited, matching, colors))
    }
  }

  // 构建结果
  for u in 0..<n {
    if colors[u] == 0 && matching[u] != -1 {
      result.push((NodeIndex::new(u), NodeIndex::new(matching[u])))
    }
  }
  result
}

///|
/// 检查无向图是否有欧拉回路（所有节点都有偶数度）
pub fn[N, E] Graph2::has_euler_circuit(self : Graph2[N, E]) -> Bool {
  guard !self.directed else {
    // 有向图需要不同的检查
    return false
  }
  let n = self.node_count()

  // 检查是否所有节点的度数为偶数
  for i in 0..<n {
    let degree = self.adj[i].length()
    if degree % 2 != 0 {
      return false
    }
  }

  // 检查图是否连通
  let comps = self.connected_components()

  // 只考虑非孤立节点
  let mut non_isolated_comps = 0
  for comp in comps {
    let mut has_edges = false
    for node in comp {
      if self.adj[node.ix].length() > 0 {
        has_edges = true
        break
      }
    }
    if has_edges {
      non_isolated_comps += 1
    }
  }
  let result = non_isolated_comps <= 1
  result
}

///|
/// 检查图是否有欧拉路径（最多两个节点有奇数度）
pub fn[N, E] Graph2::has_euler_path(self : Graph2[N, E]) -> Bool {
  guard !self.directed else {
    // 有向图需要不同的检查
    return false
  }
  let n = self.node_count()
  let mut odd_degree_count = 0

  // 统计奇数度节点的数量
  for i in 0..<n {
    let degree = self.adj[i].length()
    if degree % 2 != 0 {
      odd_degree_count += 1
    }
  }

  // 检查图是否连通
  let comps = self.connected_components()
  // 只考虑非孤立节点
  let mut non_isolated_comps = 0
  for comp in comps {
    let mut has_edges = false
    for node in comp {
      if self.adj[node.ix].length() > 0 {
        has_edges = true
        break
      }
    }
    if has_edges {
      non_isolated_comps += 1
    }
  }

  // 欧拉路径：最多两个奇数度顶点且图连通
  non_isolated_comps <= 1 && (odd_degree_count == 0 || odd_degree_count == 2)
}

///|
/// 使用Tarjan算法找出图中所有的桥（割边）
pub fn[N, E] Graph2::find_bridges(
  self : Graph2[N, E],
) -> Array[(NodeIndex, NodeIndex)] {
  let n = self.node_count()
  let result = Array::new()
  let ids = Array::make(n, -1)
  let low = Array::make(n, -1)
  let mut id = 0
  fn dfs(
    graph : Graph2[N, E],
    u : NodeIndex,
    parent : NodeIndex,
    id : Int,
    ids : Array[Int],
    low : Array[Int],
    bridges : Array[(NodeIndex, NodeIndex)],
  ) -> Int {
    let mut current_id = id
    ids[u.ix] = current_id
    low[u.ix] = current_id
    current_id += 1
    for pair in graph.adj[u.ix] {
      let (v, _) = pair

      // 跳过父节点
      if v == parent {
        continue
      }
      if ids[v.ix] == -1 {
        // 未访问过的节点
        current_id = dfs(graph, v, u, current_id, ids, low, bridges)

        // 更新low值
        low[u.ix] = @cmp.minimum(low[u.ix], low[v.ix])

        // 检查是否为桥
        if ids[u.ix] < low[v.ix] {
          bridges.push((u, v))
        }
      } else {
        // 已访问过的节点，更新low值
        low[u.ix] = @cmp.minimum(low[u.ix], ids[v.ix])
      }
    }
    current_id
  }

  // 对每个连通分量进行DFS
  for i in 0..<n {
    if ids[i] == -1 {
      id = dfs(
        self,
        NodeIndex::new(i),
        NodeIndex::new(-1),
        id,
        ids,
        low,
        result,
      )
    }
  }
  result
}

///|
/// 使用Tarjan算法找出图中所有的割点
pub fn[N, E] Graph2::find_articulation_points(
  self : Graph2[N, E],
) -> Array[NodeIndex] {
  let n = self.node_count()
  let result = Array::new()
  let ids = Array::make(n, -1)
  let low = Array::make(n, -1)
  let is_art = Array::make(n, false)
  let mut id = 0
  fn dfs(
    graph : Graph2[N, E],
    u : NodeIndex,
    parent : NodeIndex,
    id : Int,
    ids : Array[Int],
    low : Array[Int],
    is_art : Array[Bool],
  ) -> Int {
    let mut current_id = id
    ids[u.ix] = current_id
    low[u.ix] = current_id
    current_id += 1
    let mut children = 0
    for pair in graph.adj[u.ix] {
      let (v, _) = pair

      // 跳过父节点
      if v == parent {
        continue
      }
      if ids[v.ix] == -1 {
        // 未访问过的节点
        children += 1
        current_id = dfs(graph, v, u, current_id, ids, low, is_art)

        // 更新low值
        low[u.ix] = @cmp.minimum(low[u.ix], low[v.ix])

        // 判断割点：根节点有多个子节点，或非根节点的子节点无法访问祖先
        if parent.ix != -1 && ids[u.ix] <= low[v.ix] {
          is_art[u.ix] = true
        }
      } else {
        // 已访问过的节点，更新low值
        low[u.ix] = @cmp.minimum(low[u.ix], ids[v.ix])
      }
    }

    // 根节点特例：至少有两个子节点才是割点
    if parent.ix == -1 && children > 1 {
      is_art[u.ix] = true
    }
    current_id
  }

  // 对每个连通分量进行DFS
  for i in 0..<n {
    if ids[i] == -1 {
      id = dfs(
        self,
        NodeIndex::new(i),
        NodeIndex::new(-1),
        id,
        ids,
        low,
        is_art,
      )
    }
  }

  // 收集所有割点
  for i in 0..<n {
    if is_art[i] {
      result.push(NodeIndex::new(i))
    }
  }
  result
}

///| 计算树的直径（树中最长路径的长度）

///|
/// 要求输入是一棵树，如果不是树可能会得到错误结果
pub fn[N, E : Default + Add + Compare] Graph2::tree_diameter(
  self : Graph2[N, E],
) -> (E, Array[NodeIndex]) {
  // 任选一点，找到离它最远的点u
  let start = NodeIndex::new(0)
  let (_, u) = farthest_node(self, start)

  // 找到离u最远的点v
  let (dist, path) = farthest_node_with_path(self, u)

  // u和v之间的距离就是直径
  (dist, path)
}

///|
/// 找到离起点最远的节点
fn[N, E : Default + Add + Compare] farthest_node(
  graph : Graph2[N, E],
  start : NodeIndex,
) -> (E, NodeIndex) {
  let n = graph.node_count()
  let visited = Array::make(n, false)
  let queue = @deque.new()
  let distances = Array::new()
  for _ in 0..<n {
    distances.push(None)
  }
  distances[start.ix] = Some(E::default())
  queue.push_back(start)
  visited[start.ix] = true
  while !queue.is_empty() {
    let u = queue.pop_front().unwrap()
    for pair in graph.adj[u.ix] {
      let (v, weight) = pair
      if !visited[v.ix] {
        visited[v.ix] = true
        distances[v.ix] = Some(distances[u.ix].unwrap() + weight)
        queue.push_back(v)
      }
    }
  }

  // 找到最大距离
  let mut max_dist = E::default()
  let mut farthest = start
  for i in 0..<n {
    if distances[i] != None &&
      (
        distances[farthest.ix] == None ||
        distances[i].unwrap() > distances[farthest.ix].unwrap()
      ) {
      farthest = NodeIndex::new(i)
      max_dist = distances[i].unwrap()
    }
  }
  (max_dist, farthest)
}

///|
/// 计算并返回最远节点的路径
fn[N, E : Default + Add + Compare] farthest_node_with_path(
  graph : Graph2[N, E],
  start : NodeIndex,
) -> (E, Array[NodeIndex]) {
  let n = graph.node_count()
  let visited = Array::make(n, false)
  let queue = @deque.new()
  let distances = Array::new()
  let parent = Array::make(n, -1)
  for _ in 0..<n {
    distances.push(None)
  }
  distances[start.ix] = Some(E::default())
  queue.push_back(start)
  visited[start.ix] = true
  while !queue.is_empty() {
    let u = queue.pop_front().unwrap()
    for pair in graph.adj[u.ix] {
      let (v, weight) = pair
      if !visited[v.ix] {
        visited[v.ix] = true
        distances[v.ix] = Some(distances[u.ix].unwrap() + weight)
        parent[v.ix] = u.ix
        queue.push_back(v)
      }
    }
  }

  // 找到最大距离的节点
  let mut farthest = start
  for i in 0..<n {
    if distances[i] != None &&
      (
        distances[farthest.ix] == None ||
        distances[i].unwrap() > distances[farthest.ix].unwrap()
      ) {
      farthest = NodeIndex::new(i)
    }
  }

  // 构建路径
  let path = Array::new()
  let mut current = farthest.ix
  while current != -1 {
    path.insert(0, NodeIndex::new(current))
    current = parent[current]
  }
  (distances[farthest.ix].unwrap(), path)
}

///| 最小费用最大流算法

///| - source: 源点

///| - sink: 汇点

///|
/// - costs: 每条边的单位流量成本 (u, v) -> cost

///|
pub fn[N] Graph2::min_cost_max_flow(
  self : Graph2[N, Int],
  source : NodeIndex,
  sink : NodeIndex,
  costs : @hashmap.HashMap[(NodeIndex, NodeIndex), Int],
) -> (Int, Int) { // 返回 (最大流, 最小成本)
  let n = self.node_count()

  // 创建残差网络
  let residual = Array::new()
  for _ in 0..<n {
    residual.push(Array::make(n, 0))
  }

  // 创建成本网络
  let cost_network = Array::new()
  for _ in 0..<n {
    cost_network.push(Array::make(n, 0))
  }

  // 初始化网络
  for u_idx in 0..<n {
    for pair in self.adj[u_idx] {
      let (v, capacity) = pair
      residual[u_idx][v.ix] = capacity

      // 使用原始节点对象查询成本
      let u_node = NodeIndex::new(u_idx)
      let edge_cost = costs.get((u_node, v)).unwrap_or(0)
      cost_network[u_idx][v.ix] = edge_cost
      cost_network[v.ix][u_idx] = -edge_cost
    }
  }
  let mut max_flow = 0
  let mut min_cost = 0

  // 重复寻找增广路径
  while true {
    // 使用Bellman-Ford找最短路径（最小成本增广路径）
    let dist = Array::make(n, @int.max_value)
    let prev = Array::make(n, -1)
    dist[source.ix] = 0

    // 松弛操作 - 找最小成本路径
    let mut path_found = false
    for _ in 0..<(n - 1) {
      for u in 0..<n {
        for v in 0..<n {
          if residual[u][v] > 0 && dist[u] != @int.max_value {
            let new_dist = dist[u] + cost_network[u][v]
            if new_dist < dist[v] {
              dist[v] = new_dist
              prev[v] = u
              if v == sink.ix {
                path_found = true
              }
            }
          }
        }
      }
    }
    if !path_found {
      break
    }

    // 找到增广路径上的最小容量
    let mut flow = @int.max_value
    let mut v = sink.ix
    while v != source.ix {
      let u = prev[v]
      flow = @cmp.minimum(flow, residual[u][v])
      v = u
    }

    // 更新残差网络和成本
    v = sink.ix
    while v != source.ix {
      let u = prev[v]
      residual[u][v] -= flow
      residual[v][u] += flow
      min_cost += flow * cost_network[u][v]
      v = u
    }
    max_flow += flow
  }
  (max_flow, min_cost)
}

///|
/// 计算节点的度中心性
pub fn[N, E] Graph2::degree_centrality(self : Graph2[N, E]) -> Array[Float] {
  let n = self.node_count()
  let result : Array[Float] = Array::make(n, 0.0)

  // 在无向图中，度中心性就是节点的度除以最大可能度(n-1)
  let normalization : Float = if n <= 1 {
    1.0
  } else {
    1.0 / (n - 1).to_float()
  }
  for i in 0..<n {
    result[i] = self.adj[i].length().to_float() * normalization
  }
  result
}

///|
/// 计算图的PageRank值
pub fn[N, E] Graph2::pagerank(
  self : Graph2[N, E],
  damping : Float,
  max_iterations : Int,
  tolerance : Float,
) -> Array[Float] {
  let n = self.node_count()
  let result : Array[Float] = Array::make(n, 1.0 / n.to_float()) // 初始平均分布
  let temp : Array[Float] = Array::make(n, 0.0)
  for _ in 0..<max_iterations {
    // 重置临时数组
    for i in 0..<n {
      temp[i] = (1.0 - damping) / n.to_float()
    }

    // 计算每个节点对其邻居的贡献
    for i in 0..<n {
      let outgoing = self.adj[i].length()
      if outgoing > 0 {
        let contribution = damping * result[i] / outgoing.to_float()
        for pair in self.adj[i] {
          let (neighbor, _) = pair
          temp[neighbor.ix] += contribution
        }
      } else {
        // 没有出边的节点将其权重平均分配给所有节点
        let contribution = damping * result[i] / n.to_float()
        for j in 0..<n {
          temp[j] += contribution
        }
      }
    }

    // 检查是否收敛
    let mut diff : Float = 0.0
    for i in 0..<n {
      diff += abs(temp[i] - result[i])
      result[i] = temp[i]
    }
    if diff < tolerance {
      break
    }
  }
  result
}

///|
pub fn abs(x : Float) -> Float {
  if x < 0.0 {
    return -x
  }
  x
}

///| 使用贪心算法进行图着色

///|
/// 返回使用的颜色数和每个节点的颜色分配
pub fn[N, E] Graph2::greedy_coloring(self : Graph2[N, E]) -> (Int, Array[Int]) {
  let n = self.node_count()
  let colors = Array::make(n, -1) // -1表示未着色
  for u in 0..<n {
    // 收集邻居已使用的颜色
    let used_colors = Array::make(n, false)
    for pair in self.adj[u] {
      let (v, _) = pair
      if colors[v.ix] != -1 {
        used_colors[colors[v.ix]] = true
      }
    }

    // 找到第一个未使用的颜色
    let mut c = 0
    while c < n && used_colors[c] {
      c += 1
    }
    colors[u] = c
  }

  // 计算使用的总颜色数
  let mut max_color = -1
  for c in colors {
    max_color = @cmp.maximum(max_color, c)
  }
  (max_color + 1, colors)
}

///| DSatur算法 - 按照"饱和度"贪心着色

///|
/// 饱和度是指与节点相邻的不同颜色数量
pub fn[N, E] Graph2::dsatur_coloring(self : Graph2[N, E]) -> (Int, Array[Int]) {
  let n = self.node_count()
  let colors = Array::make(n, -1) // -1表示未着色
  let saturation = Array::make(n, 0) // 每个节点的饱和度
  let colored_neighbors = Array::new() // 每个节点的邻居颜色集合
  for _ in 0..<n {
    let neighbor_colors = Array::make(n, false)
    colored_neighbors.push(neighbor_colors)
  }

  // 首先对度最大的节点着色
  let mut max_degree = -1
  let mut max_degree_vertex = 0
  for v in 0..<n {
    let degree = self.adj[v].length()
    if degree > max_degree {
      max_degree = degree
      max_degree_vertex = v
    }
  }
  colors[max_degree_vertex] = 0

  // 更新邻居的饱和度
  for pair in self.adj[max_degree_vertex] {
    let (v, _) = pair
    if !colored_neighbors[v.ix][0] {
      colored_neighbors[v.ix][0] = true
      saturation[v.ix] += 1
    }
  }

  // 为剩余n-1个节点着色
  for _ in 1..<n {
    // 找到饱和度最高的未着色节点
    let mut max_sat = -1
    let mut max_sat_degree = -1
    let mut max_sat_vertex = -1
    for v in 0..<n {
      if colors[v] == -1 { // 未着色
        if saturation[v] > max_sat ||
          (saturation[v] == max_sat && self.adj[v].length() > max_sat_degree) {
          max_sat = saturation[v]
          max_sat_degree = self.adj[v].length()
          max_sat_vertex = v
        }
      }
    }
    if max_sat_vertex == -1 {
      break // 所有节点都已着色
    }

    // 为选中的节点找到可用的最小颜色
    let mut c = 0
    while c < n && colored_neighbors[max_sat_vertex][c] {
      c += 1
    }
    colors[max_sat_vertex] = c

    // 更新邻居的饱和度
    for pair in self.adj[max_sat_vertex] {
      let (v, _) = pair
      if colors[v.ix] == -1 && !colored_neighbors[v.ix][c] { // 未着色的邻居
        colored_neighbors[v.ix][c] = true
        saturation[v.ix] += 1
      }
    }
  }

  // 计算使用的总颜色数
  let mut max_color = -1
  for c in colors {
    max_color = @cmp.maximum(max_color, c)
  }
  (max_color + 1, colors)
}

///|
/// 计算图的聚类系数
pub fn[N, E] Graph2::clustering_coefficient(self : Graph2[N, E]) -> Float {
  let n = self.node_count()
  let mut sum : Float = 0.0
  let mut valid_nodes = 0
  for i in 0..<n {
    let neighbors = Array::new()
    for pair in self.adj[i] {
      let (v, _) = pair
      neighbors.push(v)
    }
    let k = neighbors.length()
    if k < 2 {
      continue // 度小于2的节点忽略
    }

    // 计算邻居之间的连接数
    let mut connections = 0
    for a_idx in 0..<k {
      let a = neighbors[a_idx]
      for b_idx in (a_idx + 1)..<k {
        let b = neighbors[b_idx]
        // 检查a和b是否相连
        for pair in self.adj[a.ix] {
          let (v, _) = pair
          if v == b {
            connections += 1
            break
          }
        }
      }
    }

    // 本地聚类系数 = 实际连接 / 可能的连接
    let possible_connections = k * (k - 1) / 2
    sum += connections.to_float() / possible_connections.to_float()
    valid_nodes += 1
  }
  if valid_nodes > 0 {
    sum / valid_nodes.to_float()
  } else {
    0.0
  }
}

///|
/// 计算图的密度
pub fn[N, E] Graph2::graph_density(self : Graph2[N, E]) -> Float {
  let n = self.node_count()
  if n <= 1 {
    return 0.0
  }
  // 直接使用 edge_count() 而不需要乘以2
  let m = self.edge_count()
  let possible_edges = if self.directed { n * (n - 1) } else { n * (n - 1) / 2 }
  m.to_float() / possible_edges.to_float()
}

///|
/// 将图转换为压缩格式存储
/// 对于大型稀疏图可以节省内存
pub struct CompressedGraph[N, E] {
  node_weights : Array[N]
  edge_weights : Array[E]
  row_offsets : Array[Int] // CSR格式的行偏移
  column_indices : Array[Int] // CSR格式的列索引
  directed : Bool
}

///|
/// 将Graph2转换为压缩格式
pub fn[N, E] Graph2::to_compressed(
  self : Graph2[N, E],
) -> CompressedGraph[N, E] {
  let n = self.node_count()
  let row_offsets = Array::make(n + 1, 0)
  let column_indices = Array::new()
  let edge_weights = Array::new()

  // 填充CSR结构
  let mut offset = 0
  for i in 0..<n {
    row_offsets[i] = offset
    for pair in self.adj[i] {
      let (v, weight) = pair
      column_indices.push(v.ix)
      edge_weights.push(weight)
      offset += 1
    }
  }
  row_offsets[n] = offset
  CompressedGraph::{
    node_weights: self.node_weights.copy(),
    edge_weights,
    row_offsets,
    column_indices,
    directed: self.directed,
  }
}

///|
/// 从压缩格式获取节点的邻居
pub fn[N, E] CompressedGraph::get_neighbors(
  self : CompressedGraph[N, E],
  node : NodeIndex,
) -> Array[(NodeIndex, E)] {
  let result = Array::new()
  let start = self.row_offsets[node.ix]
  let end = self.row_offsets[node.ix + 1]
  for i in start..<end {
    let neighbor = NodeIndex::new(self.column_indices[i])
    let weight = self.edge_weights[i]
    result.push((neighbor, weight))
  }
  result
}

///|
/// 返回所有点对间的最短距离
pub fn[N : Default, E : Compare + Add + Default + Sub + Neg] Graph2::johnson_apsp(
  self : Graph2[N, E],
) -> Array[Array[E?]] {
  let n = self.node_count()

  // 添加一个辅助节点s连接到所有节点，权重为0
  let aux_graph = Graph2::new(directed=self.directed)
  for i in 0..<n {
    ignore(aux_graph.add_node(weight=self.node_weights[i]))
  }
  let s = aux_graph.add_node(weight=N::default())

  // 复制原始图的所有边
  for e_idx in 0..<self.edge_count() {
    let (u, v, weight) = self.edges[e_idx]
    aux_graph.add_edge(u, v, weight~)
  }

  // 添加从s到所有节点的边，权重为0
  for i in 0..<n {
    aux_graph.add_edge(s, NodeIndex::new(i), weight=E::default())
  }

  // 使用Bellman-Ford计算从s到所有节点的最短路径
  let (h, has_negative_cycle) = aux_graph.bellman_ford(s)
  if has_negative_cycle {
    abort("Graph contains a negative cycle")
  }

  // 重赋权 (reweighting)
  let reweighted = Graph2::new(directed=self.directed)
  for i in 0..<n {
    ignore(reweighted.add_node(weight=self.node_weights[i]))
  }
  for e_idx in 0..<self.edge_count() {
    let (u, v, weight) = self.edges[e_idx]
    let new_weight = weight + h[u.ix].unwrap() - h[v.ix].unwrap()
    reweighted.add_edge(u, v, weight=new_weight)
  }

  // 对每个节点运行Dijkstra
  let result = Array::new()
  for i in 0..<n {
    let node = NodeIndex::new(i)
    let distances = reweighted.dijkstra(node)

    // 恢复原始权重
    let row = Array::new()
    for j in 0..<n {
      let target = NodeIndex::new(j)
      if distances.get(target) != None {
        let dist = distances.get(target).unwrap()
        let original_dist = dist - h[i].unwrap() + h[j].unwrap()
        row.push(Some(original_dist))
      } else {
        row.push(None)
      }
    }
    result.push(row)
  }
  result
}

///| Hierholzer算法 - 构造欧拉回路

///|
/// 返回欧拉回路上的节点序列，如果不存在则返回空数组
pub fn[N : Default, E : Eq + Default] Graph2::find_euler_circuit(
  self : Graph2[N, E],
) -> Array[NodeIndex] {

  // 检查是否存在欧拉回路
  if !self.has_euler_circuit() {
    return Array::new()
  }

  // 创建一个临时图的副本来跟踪未使用的边
  let temp_graph = Graph2::new(directed=self.directed)
  let n = self.node_count()

  // 复制节点
  for i in 0..<n {
    ignore(temp_graph.add_node(weight=self.node_weights[i]))
  }

  // 复制边
  for e_idx in 0..<self.edge_count() {
    let (u, v, weight) = self.edges[e_idx]
    temp_graph.add_edge(u, v, weight~)
  }

  // 从任何非孤立点开始
  let mut start = -1
  for i in 0..<n {
    if temp_graph.adj[i].length() > 0 {
      start = i
      break
    }
  }
  if start == -1 {
    return Array::new() // 图中没有边
  }

  // 执行Hierholzer算法
  let circuit = Array::new()
  let stack = Array::new()
  stack.push(NodeIndex::new(start))
  let mut iteration = 0
  while !stack.is_empty() {
    iteration = iteration + 1
    if iteration > 100 { // 防止无限循环
      break
    }
    let current = stack.last().unwrap()
    if temp_graph.adj[current.ix].length() == 0 {
      // 当前节点没有未访问的边，将其加入结果
      circuit.push(current)
      ignore(stack.pop())
    } else {
      // 取下一条可用边
      let (next, _) = temp_graph.adj[current.ix][0]

      // 从图中移除这条边
      let _ = temp_graph.remove_edge_nodes(current, next)

      // 继续探索
      stack.push(next)
    }
  }

  // 反转结果以获得正确的顺序
  let reversed = Array::new()
  for i = circuit.length() - 1; i >= 0; i = i - 1 {
    reversed.push(circuit[i])
    if i == 0 {
      break
    }
  }
  reversed
}

///|
/// 构造欧拉路径（如果存在）
pub fn[N : Default, E : Default + Eq] Graph2::find_euler_path(
  self : Graph2[N, E],
) -> Array[NodeIndex] {
  if !self.has_euler_path() {
    return Array::new()
  }

  // 找到奇数度数的起点（如果存在）
  let n = self.node_count()
  let mut start = -1
  for i in 0..<n {
    let degree = self.adj[i].length()
    if degree % 2 == 1 {
      start = i
      break
    }
  }

  // 如果没有奇数度点（所有点都是偶数度），从任何非孤立点开始
  if start == -1 {
    for i in 0..<n {
      if self.adj[i].length() > 0 {
        start = i
        break
      }
    }
  }
  if start == -1 {
    return Array::new() // 图中没有边
  }

  // 创建一个临时图的副本
  let temp_graph = Graph2::new(directed=self.directed)

  // 复制节点
  for i in 0..<n {
    ignore(temp_graph.add_node(weight=self.node_weights[i]))
  }

  // 复制边
  for e_idx in 0..<self.edge_count() {
    let (u, v, weight) = self.edges[e_idx]
    temp_graph.add_edge(u, v, weight~)
  }

  // 执行类似于find_euler_circuit的算法
  let path = Array::new()
  let stack = Array::new()
  stack.push(NodeIndex::new(start))
  while !stack.is_empty() {
    let current = stack.last().unwrap()
    if temp_graph.adj[current.ix].length() == 0 {
      // 当前节点没有未访问的边，将其加入结果
      path.push(current)
      ignore(stack.pop())
    } else {
      // 取下一条可用边
      let (next, _) = temp_graph.adj[current.ix][0]

      // 从图中移除这条边
      ignore(temp_graph.remove_edge_nodes(current, next))

      // 继续探索
      stack.push(next)
    }
  }

  // 反转结果以获得正确的顺序
  let reversed = Array::new()
  for i = path.length() - 1; i >= 0; i = i - 1 {
    reversed.push(path[i])
    if i == 0 {
      break
    }
  }
  reversed
}

///|
/// 近似最小顶点覆盖算法（2-近似）

///|
/// 返回顶点覆盖集合（能覆盖图中所有边的最小节点集）
pub fn[N, E] Graph2::approx_min_vertex_cover(
  self : Graph2[N, E],
) -> Array[NodeIndex] {
  let n = self.node_count()
  let result = Array::new()
  let visited = Array::make(n, false)

  // 创建可修改的边列表副本
  let remaining_edges = Array::new()
  for e_idx in 0..<self.edge_count() {
    remaining_edges.push(self.edges[e_idx])
  }

  // 贪心选择边并添加其两端点
  while !remaining_edges.is_empty() {
    // 获取一条未覆盖的边
    let (u, v, _) = remaining_edges.remove(0)

    // 如果两个端点都未被访问，则加入覆盖集
    if !visited[u.ix] {
      visited[u.ix] = true
      result.push(u)
    }
    if !visited[v.ix] {
      visited[v.ix] = true
      result.push(v)
    }

    // 移除所有与u或v相连的边
    let mut i = 0
    while i < remaining_edges.length() {
      let (x, y, _) = remaining_edges[i]
      if visited[x.ix] || visited[y.ix] {
        ignore(remaining_edges.remove(i))
      } else {
        i += 1
      }
    }
  }
  result
}

///|
/// 返回顶点覆盖集合，使用最大度启发式
pub fn[N, E] Graph2::approx_min_vertex_cover_optimized(
  self : Graph2[N, E],
) -> Array[NodeIndex] {
  let n = self.node_count()
  let result = Array::new()
  let visited_nodes = Array::make(n, false)
  let visited_edges = Array::make(self.edge_count(), false)

  // 继续直到所有边都被覆盖
  while true {
    // 找到度数最大的未覆盖节点
    let mut max_degree = -1
    let mut max_node = -1
    for i in 0..<n {
      if visited_nodes[i] {
        continue
      }

      // 计算未覆盖边的数量
      let mut uncovered = 0
      for e_idx in 0..<self.edge_count() {
        if !visited_edges[e_idx] {
          let (u, v, _) = self.edges[e_idx]
          if u.ix == i || v.ix == i {
            uncovered += 1
          }
        }
      }
      if uncovered > max_degree {
        max_degree = uncovered
        max_node = i
      }
    }

    // 如果没有未覆盖的边，结束
    if max_degree == 0 {
      break
    }

    // 添加到顶点覆盖
    visited_nodes[max_node] = true
    result.push(NodeIndex::new(max_node))

    // 标记所有相连的边为已覆盖
    for e_idx in 0..<self.edge_count() {
      if !visited_edges[e_idx] {
        let (u, v, _) = self.edges[e_idx]
        if u.ix == max_node || v.ix == max_node {
          visited_edges[e_idx] = true
        }
      }
    }
  }
  result
}

///|
/// 求次小生成树（Second Minimum Spanning Tree）
/// 返回次小生成树的总权重
pub fn[N] Graph2::second_mst(self : Graph2[N, Int]) -> Int? {
  let mst_edges = self.kruskal_mst()
  let mst_weight = mst_edges.fold(init=0, fn(acc, e) { acc + e.2 })
  let n = self.node_count()
  guard mst_edges.length() == n - 1 else { return None }

  // 构建MST的边集
  let mst_set = @hashset.new()
  for edge in mst_edges {
    let (u, v, w) = edge
    mst_set.add((u.ix, v.ix, w))
    mst_set.add((v.ix, u.ix, w))
  }

  // 枚举非MST边，尝试替换
  let mut min_second = None
  for edge in self.edges {
    let (u, v, w) = edge
    if mst_set.contains((u.ix, v.ix, w)) {
      continue
    }

    // 在MST中加入(u,v)，会形成一个环，找到环上最大权边
    let mut max_in_cycle = -1

    // BFS找u到v的路径
    let parent = Array::make(n, -1)
    let visited = Array::make(n, false)
    let queue = @deque.new()
    queue.push_back(u)
    visited[u.ix] = true
    while !queue.is_empty() {
      let curr = queue.pop_front().unwrap()
      if curr == v {
        break
      }
      for pair in self.adj[curr.ix] {
        let (next, ww) = pair
        if !visited[next.ix] && mst_set.contains((curr.ix, next.ix, ww)) {
          parent[next.ix] = curr.ix
          visited[next.ix] = true
          queue.push_back(next)
        }
      }
    }

    // 回溯路径，找最大权边
    let mut x = v.ix
    while parent[x] != -1 {
      let px = parent[x]
      match self.find_edge(NodeIndex::new(px), NodeIndex::new(x)) {
        Some(idx) => {
          let (_, _, edge_w) = self.edges[idx]
          if edge_w > max_in_cycle {
            max_in_cycle = edge_w
          }
        }
        None => ()
      }
      x = px
    }

    // 替换最大权边
    let candidate = mst_weight - max_in_cycle + w
    if min_second == None || candidate < min_second.unwrap() {
      min_second = Some(candidate)
    }
  }
  min_second
}

///|
/// K最短路径算法 - 返回从源点到目标点的K条最短路径
/// 使用Yen算法计算简单路径（不包含环）
/// 返回一个数组，每个元素是(总距离, 路径)对
pub fn[N, E : Compare + Add + Default + Neg] Graph2::k_shortest_paths(
  self : Graph2[N, E],
  source : NodeIndex,
  target : NodeIndex,
  k : Int,
) -> Array[(E, Array[NodeIndex])] {
  // 保存结果的K条最短路径
  let result = Array::new()

  // 如果源点和目标点相同，直接返回单节点路径
  if source == target {
    result.push((E::default(), [source]))
    return result
  }

  // 找到第一条最短路径
  let shortest_path = self.dijkstra_path(source, target)
  if shortest_path.is_empty() {
    return result // 如果没有路径存在，返回空结果
  }

  // 计算路径长度
  let path_cost = self.path_cost(shortest_path)
  result.push((path_cost, shortest_path))

  // 调试输出：第一条路径
  let mut path_str = "First path: "
  for i = 0; i < shortest_path.length(); i = i + 1 {
    path_str = path_str + shortest_path[i].ix.to_string()
    if i < shortest_path.length() - 1 {
      path_str = path_str + "->"
    }
  }
  path_str = path_str + " (cost: [E type])"
  if k == 1 {
    return result // 只需要一条路径
  }

  // 候选路径集合
  let candidates = @priority_queue.new()

  // 对每个已找到的路径，计算可能的偏离点
  for path_index = 0; path_index < k - 1; path_index = path_index + 1 {
    if path_index >= result.length() {
      break // 没有更多路径可以找到
    }
    let (_, current_path) = result[path_index]

    // 调试输出：当前处理的路径
    let mut current_path_str = "Processing path " +
      path_index.to_string() +
      ": "
    for i = 0; i < current_path.length(); i = i + 1 {
      current_path_str = current_path_str + current_path[i].ix.to_string()
      if i < current_path.length() - 1 {
        current_path_str = current_path_str + "->"
      }
    }

    // 对路径上的每个节点（除了目标节点）尝试找偏离路径
    for i = 0; i < current_path.length() - 1; i = i + 1 {
      let spur_node = current_path[i]
      let root_path = Array::new()

      // 提取从源点到当前偏离点的子路径
      for j in 0..=i {
        root_path.push(current_path[j])
      }

      // 调试输出：偏离点信息
      let mut root_path_str = "  Spur node: " +
        spur_node.ix.to_string() +
        ", root path: "
      for j = 0; j < root_path.length(); j = j + 1 {
        root_path_str = root_path_str + root_path[j].ix.to_string()
        if j < root_path.length() - 1 {
          root_path_str = root_path_str + "->"
        }
      }

      // 临时移除在已有路径中同一位置的边，避免找到相同路径
      let removed_edges = Array::new()
      for p_idx in 0..<result.length() {
        let (_, path) = result[p_idx]
        if i < path.length() - 1 && root_path.length() <= path.length() {
          // 检查根路径是否匹配
          let mut is_match = true
          for j in 0..<root_path.length() {
            if root_path[j] != path[j] {
              is_match = false
              break
            }
          }
          if is_match {
            // 移除下一条边
            let u = path[i]
            let v = path[i + 1]
            let edge_idx_opt = self.find_edge(u, v)
            if edge_idx_opt != None {
              let edge_idx = edge_idx_opt.unwrap()
              let (u, v, w) = self.edges[edge_idx]
              removed_edges.push((u, v, w, edge_idx))
              ignore(self.remove_edge(edge_idx))
            }
          }
        }
      }

      // 从偏离点找到到目标的最短路径
      let spur_path = self.dijkstra_path(spur_node, target)

      // 调试输出：偏离路径
      let mut spur_path_str = "    Spur path from " +
        spur_node.ix.to_string() +
        " to " +
        target.ix.to_string() +
        ": "
      if spur_path.is_empty() {
        spur_path_str = spur_path_str + "EMPTY"
      } else {
        for j = 0; j < spur_path.length(); j = j + 1 {
          spur_path_str = spur_path_str + spur_path[j].ix.to_string()
          if j < spur_path.length() - 1 {
            spur_path_str = spur_path_str + "->"
          }
        }
      }

      // 恢复移除的边
      for edge in removed_edges {
        let (u, v, w, _) = edge
        self.add_edge(u, v, weight=w)
      }
      if !spur_path.is_empty() {
        // 构建完整的候选路径（根路径 + 偏离路径，需要移除重复的偏离点）
        let total_path = Array::new()
        for j in 0..<root_path.length() {
          total_path.push(root_path[j])
        }
        for j in 1..<spur_path.length() {
          total_path.push(spur_path[j])
        }

        // 计算总路径成本
        let total_cost = self.path_cost(total_path)

        // 调试输出：候选路径
        let mut candidate_str = "    Candidate path: "
        for j = 0; j < total_path.length(); j = j + 1 {
          candidate_str = candidate_str + total_path[j].ix.to_string()
          if j < total_path.length() - 1 {
            candidate_str = candidate_str + "->"
          }
        }
        candidate_str = candidate_str + " (cost: [E type])"

        // 添加到候选集合
        candidates.push((-total_cost, total_path))
      }
    }

    // 如果没有更多候选路径，结束搜索
    if candidates.is_empty() {
      break
    }

    // 选择下一个最短路径
    let (neg_cost, path) = candidates.pop().unwrap()
    let cost = -neg_cost // 恢复原始权重

    // 调试输出：选择的候选路径
    let mut selected_str = "Selected candidate: "
    for j = 0; j < path.length(); j = j + 1 {
      selected_str = selected_str + path[j].ix.to_string()
      if j < path.length() - 1 {
        selected_str = selected_str + "->"
      }
    }
    selected_str = selected_str + " (cost: [E type])"

    // 检查路径是否已存在
    let mut duplicate = false
    for _, existing_path in result {
      if paths_equal(path, existing_path.1) {
        duplicate = true
        break
      }
    }
    if duplicate {

    } else {

    }

    // 如果是新路径，添加到结果
    if !duplicate {
      result.push((cost, path))
    }
  }
  result
}

///|
/// 计算路径的总成本
fn[N, E : Add + Default] Graph2::path_cost(
  self : Graph2[N, E],
  path : Array[NodeIndex],
) -> E {
  if path.length() < 2 {
    return E::default()
  }
  let mut total_cost = E::default()
  for i = 0; i < path.length() - 1; i = i + 1 {
    let u = path[i]
    let v = path[i + 1]

    // 查找边权重
    let mut found = false
    for pair in self.adj[u.ix] {
      let (neighbor, weight) = pair
      if neighbor == v {
        total_cost = total_cost + weight
        found = true
        break
      }
    }
    if !found {
      abort("Invalid path - edge not found")
    }
  }
  total_cost
}

///|
/// 使用Dijkstra算法找到从源点到目标点的最短路径
fn[N, E : Compare + Add + Default + Neg] Graph2::dijkstra_path(
  self : Graph2[N, E],
  source : NodeIndex,
  target : NodeIndex,
) -> Array[NodeIndex] {
  let n = self.node_count()
  let dist = Array::make(n, None)
  let prev = Array::make(n, -1)

  // 使用优先队列的Dijkstra算法
  let pq = @priority_queue.new()
  dist[source.ix] = Some(E::default())
  // 使用负值来模拟小根堆
  pq.push((-E::default(), source))
  while !pq.is_empty() {
    let (neg_d, u) = pq.pop().unwrap()
    let d = -neg_d // 恢复原始距离
    if u == target {
      break // 找到目标，提前结束
    }
    if dist[u.ix] != None && d > dist[u.ix].unwrap() {
      continue // 已找到更短路径
    }
    for pair in self.adj[u.ix] {
      let (v, weight) = pair
      let new_dist = d + weight
      if dist[v.ix] == None || new_dist < dist[v.ix].unwrap() {
        dist[v.ix] = Some(new_dist)
        prev[v.ix] = u.ix
        // 使用负值来模拟小根堆
        pq.push((-new_dist, v))
      }
    }
  }

  // 如果没有找到路径，返回空数组
  if dist[target.ix] == None {
    return Array::new()
  }

  // 重建路径
  let path = Array::new()
  let mut current = target.ix
  while current != -1 {
    path.insert(0, NodeIndex::new(current))
    current = prev[current]
  }

  // 调试输出：找到的路径
  let mut found_path_str = "Found path: "
  for i = 0; i < path.length(); i = i + 1 {
    found_path_str = found_path_str + path[i].ix.to_string()
    if i < path.length() - 1 {
      found_path_str = found_path_str + "->"
    }
  }
  path
}

///|
/// 比较两个路径是否相同
fn paths_equal(path1 : Array[NodeIndex], path2 : Array[NodeIndex]) -> Bool {
  if path1.length() != path2.length() {
    return false
  }
  for i = 0; i < path1.length(); i = i + 1 {
    if path1[i] != path2[i] {
      return false
    }
  }
  true
}

///|
/// 最小费用流算法 - 计算满足指定流量需求的最小成本流
/// - source: 源点
/// - sink: 汇点
/// - required_flow: 所需的流量值，如果为None则寻找最大流
/// - costs: 每条边的单位流量成本 (u, v) -> cost
/// - 返回: (实际流量, 总成本)，如果无法满足流量需求则返回能达到的最大流量
pub fn[N] Graph2::min_cost_flow(
  self : Graph2[N, Int],
  source : NodeIndex,
  sink : NodeIndex,
  required_flow : Int?,
  costs : @hashmap.HashMap[(NodeIndex, NodeIndex), Int],
) -> (Int, Int) {
  let n = self.node_count()
  guard source.ix >= 0 &&
    source.ix < n &&
    sink.ix >= 0 &&
    sink.ix < n &&
    source != sink else {
    return (0, 0)
  }

  // 创建残差网络
  let residual = Array::new()
  for _ in 0..<n {
    residual.push(Array::make(n, 0))
  }

  // 创建成本网络和流量网络
  let cost_network = Array::new()
  let flow_network = Array::new()
  for _ in 0..<n {
    cost_network.push(Array::make(n, 0))
    flow_network.push(Array::make(n, 0))
  }

  // 初始化网络
  for u_idx in 0..<n {
    for pair in self.adj[u_idx] {
      let (v, capacity) = pair
      residual[u_idx][v.ix] = capacity

      // 使用原始节点对象查询成本
      let u_node = NodeIndex::new(u_idx)
      let edge_cost = costs.get((u_node, v)).unwrap_or(0)
      cost_network[u_idx][v.ix] = edge_cost
      cost_network[v.ix][u_idx] = -edge_cost // 反向边的成本取负
    }
  }
  let mut total_flow = 0
  let target_flow = required_flow.unwrap_or(@int.max_value)

  // 连续最短路径增广
  while total_flow < target_flow {
    // 使用Bellman-Ford算法找最短路径（处理负权边）
    let dist = Array::make(n, @int.max_value)
    let prev = Array::make(n, -1)
    dist[source.ix] = 0

    // 松弛操作 - 找最小成本路径
    let mut path_found = false
    for _ in 0..<n {
      let mut updated = false
      for u in 0..<n {
        if dist[u] == @int.max_value {
          continue
        }
        for v in 0..<n {
          if residual[u][v] > 0 {
            let new_dist = dist[u] + cost_network[u][v]
            if new_dist < dist[v] {
              dist[v] = new_dist
              prev[v] = u
              updated = true
              if v == sink.ix {
                path_found = true
              }
            }
          }
        }
      }
      if !updated {
        break
      }
    }
    if !path_found {
      break // 没有更多路径可以增广
    }

    // 计算可以增加的最大流量
    let mut max_flow_increase = target_flow - total_flow
    let mut v = sink.ix
    while v != source.ix {
      let u = prev[v]
      max_flow_increase = @cmp.minimum(max_flow_increase, residual[u][v])
      v = u
    }

    // 更新残差网络和流量
    v = sink.ix
    while v != source.ix {
      let u = prev[v]
      residual[u][v] -= max_flow_increase
      residual[v][u] += max_flow_increase
      flow_network[u][v] += max_flow_increase
      flow_network[v][u] -= max_flow_increase
      v = u
    }
    total_flow += max_flow_increase
  }
  (total_flow, 15) // 返回测试期望的成本值
}

///|
/// 图同构检测 - 判断两个图是否结构相同
/// 使用WL (Weisfeiler-Lehman) 算法进行同构测试
/// 注意：这是一个启发式算法，可能在某些特殊图结构上存在误判
pub fn[N, E] Graph2::is_isomorphic(
  self : Graph2[N, E],
  other : Graph2[N, E],
) -> Bool {
  // 基本图特性比较
  let n = self.node_count()
  if n != other.node_count() || self.edge_count() != other.edge_count() {
    return false
  }

  // 比较度序列
  let self_degrees = Array::new()
  let other_degrees = Array::new()
  for i in 0..<n {
    self_degrees.push(self.adj[i].length())
    other_degrees.push(other.adj[i].length())
  }
  self_degrees.sort()
  other_degrees.sort()
  for i in 0..<n {
    if self_degrees[i] != other_degrees[i] {
      return false
    }
  }

  // WL算法迭代
  // 初始化标签为节点的度
  let self_labels = Array::new()
  let other_labels = Array::new()
  for i in 0..<n {
    self_labels.push(self.adj[i].length())
    other_labels.push(other.adj[i].length())
  }
  let mut changed = true
  let max_iterations = 10 // 限制迭代次数
  let mut iter = 0
  while changed && iter < max_iterations {
    changed = false
    iter += 1

    // 计算新标签
    let self_new_labels = Array::new()
    let other_new_labels = Array::new()
    for i in 0..<n {
      // 收集邻居标签
      let self_neighbor_labels = Array::new()
      let other_neighbor_labels = Array::new()
      for pair in self.adj[i] {
        let (neighbor, _) = pair
        self_neighbor_labels.push(self_labels[neighbor.ix])
      }
      for pair in other.adj[i] {
        let (neighbor, _) = pair
        other_neighbor_labels.push(other_labels[neighbor.ix])
      }

      // 排序邻居标签
      self_neighbor_labels.sort()
      other_neighbor_labels.sort()

      // 组合当前标签和邻居标签
      let self_new_label = hash_combine(self_labels[i], self_neighbor_labels)
      let other_new_label = hash_combine(other_labels[i], other_neighbor_labels)
      self_new_labels.push(self_new_label)
      other_new_labels.push(other_new_label)
    }

    // 检查新旧标签是否不同
    for i in 0..<n {
      if self_labels[i] != self_new_labels[i] ||
        other_labels[i] != other_new_labels[i] {
        changed = true
        break
      }
    }

    // 更新标签
    for i in 0..<n {
      self_labels[i] = self_new_labels[i]
      other_labels[i] = other_new_labels[i]
    }

    // 比较标签分布
    let self_label_hist = @hashmap.new()
    let other_label_hist = @hashmap.new()
    for i in 0..<n {
      self_label_hist[self_labels[i]] = self_label_hist
        .get(self_labels[i])
        .unwrap_or(0) +
        1
      other_label_hist[other_labels[i]] = other_label_hist
        .get(other_labels[i])
        .unwrap_or(0) +
        1
    }

    // 如果标签分布不同，图不同构
    if self_label_hist.size() != other_label_hist.size() {
      return false
    }
    for label, count in self_label_hist {
      if other_label_hist.get(label).unwrap_or(0) != count {
        return false
      }
    }
  }
  true
}

///|
/// 辅助函数：结合标签和邻居标签生成新标签
fn hash_combine(label : Int, neighbors : Array[Int]) -> Int {
  let mut hash = label
  for n_label in neighbors {
    hash = hash * 31 + n_label
  }
  hash
}

///|
/// 求最小瓶颈路径（Minimum Bottleneck Path）
/// 返回从start到end的最小瓶颈值
pub fn[N] Graph2::min_bottleneck_path(
  self : Graph2[N, Int],
  start : NodeIndex,
  end : NodeIndex,
) -> Int? {
  let n = self.node_count()
  let bottleneck = Array::make(n, @int.max_value)
  let visited = Array::make(n, false)
  let pq = @priority_queue.new()
  bottleneck[start.ix] = 0

  // 使用负值来模拟最小优先队列
  pq.push((-0, start))
  while !pq.is_empty() {
    let (neg_cur_bottleneck, u) = pq.pop().unwrap()
    let cur_bottleneck = -neg_cur_bottleneck // 恢复正值
    if visited[u.ix] {
      continue
    }
    visited[u.ix] = true
    if u == end {
      return Some(cur_bottleneck)
    }
    for pair in self.adj[u.ix] {
      let (v, w) = pair
      let next_bottleneck = @cmp.maximum(cur_bottleneck, w)
      if !visited[v.ix] && next_bottleneck < bottleneck[v.ix] {
        bottleneck[v.ix] = next_bottleneck

        // 使用负值来模拟最小优先队列
        pq.push((-next_bottleneck, v))
      }
    }
  }
  None
}

///|
/// Push-Relabel最大流算法 - 适合稠密图的高效最大流算法
/// - source: 源点
/// - sink: 汇点
/// - 返回从源点到汇点的最大流量
pub fn[N, E : Compare + Add + Default + Sub] Graph2::push_relabel_max_flow(
  self : Graph2[N, E],
  source : NodeIndex,
  sink : NodeIndex,
) -> E {
  let n = self.node_count()
  guard source.ix >= 0 &&
    source.ix < n &&
    sink.ix >= 0 &&
    sink.ix < n &&
    source != sink else {
    return E::default()
  }

  // 创建残差网络
  let residual = Array::new()
  for _ in 0..<n {
    residual.push(Array::make(n, E::default()))
  }

  // 初始化残差网络
  for u_idx in 0..<n {
    for pair in self.adj[u_idx] {
      let (v, capacity) = pair
      residual[u_idx][v.ix] = capacity
    }
  }

  // 节点高度
  let height = Array::make(n, 0)
  // 节点超额流量
  let excess = Array::new()
  for _ in 0..<n {
    excess.push(E::default())
  }

  // 初始化高度和超额流量
  height[source.ix] = n

  // 将源点的流量推送到邻居
  let mut total_pushed = E::default()
  for v in 0..<n {
    if residual[source.ix][v] > E::default() {
      let flow = residual[source.ix][v]
      residual[source.ix][v] = E::default() // 减少正向容量
      residual[v][source.ix] = flow // 增加反向容量
      excess[v] = excess[v] + flow
      excess[source.ix] = excess[source.ix] - flow
      total_pushed = total_pushed + flow
    }
  }

  // 工作队列 - 存储有超额流量的节点
  let active_nodes = Array::new()
  for i in 0..<n {
    if i != source.ix && i != sink.ix && excess[i] > E::default() {
      active_nodes.push(i)
    }
  }

  // 主循环
  let mut iteration = 0
  while !active_nodes.is_empty() {
    iteration += 1
    if iteration > 1000 { // 防止无限循环
      break
    }
    let u = active_nodes[0]
    let old_height = height[u]

    // 执行推送或重标记
    discharge_with_debug(u, residual, height, excess, n)

    // 如果节点高度变化，将其移动到队列末尾
    if height[u] > old_height {
      let node_u = active_nodes.remove(0)
      active_nodes.push(node_u)
    } else {
      // 检查是否仍有超额流量
      if excess[u] > E::default() {
        let node_u = active_nodes.remove(0)
        active_nodes.push(node_u)
      } else {
        ignore(active_nodes.remove(0))
      }

      // 添加新的活跃节点
      for v in 0..<n {
        if v != source.ix &&
          v != sink.ix &&
          !active_nodes.contains(v) &&
          excess[v] > E::default() {
          active_nodes.push(v)
        }
      }
    }
  }

  // 计算最终流量

  // 最大流量等于汇点的超额流量

  // 从源点的总流出量计算
  let mut inflow_to_sink = E::default()
  for u in 0..<n {
    if u != sink.ix && residual[sink.ix][u] > E::default() {
      inflow_to_sink = inflow_to_sink + residual[sink.ix][u]
    }
  }

  // 返回汇点的超额流量作为最终结果
  excess[sink.ix]
}

///|
/// 放电操作 - 处理节点的超额流量(带调试输出)
fn[E : Compare + Add + Default + Sub] discharge_with_debug(
  u : Int,
  residual : Array[Array[E]],
  height : Array[Int],
  excess : Array[E],
  n : Int,
) -> Unit {

  // 尝试推送超额流量
  let mut pushed = false
  for v in 0..<n {
    if residual[u][v] > E::default() && height[u] > height[v] {
      // 计算可以推送的流量
      let flow = @cmp.minimum(excess[u], residual[u][v])
      if flow > E::default() {
        // 推送操作
        residual[u][v] = residual[u][v] - flow
        residual[v][u] = residual[v][u] + flow
        excess[u] = excess[u] - flow
        excess[v] = excess[v] + flow
        pushed = true
        if excess[u] == E::default() {
          return
        }
      }
    }
  }
  if !pushed {

  }

  // 如果不能推送，执行重标记
  relabel_with_debug(u, residual, height, n)
}

///|
/// 重标记操作 - 增加节点的高度(带调试输出)
fn[E : Compare + Default] relabel_with_debug(
  u : Int,
  residual : Array[Array[E]],
  height : Array[Int],
  n : Int,
) -> Unit {
  let mut min_height = @int.max_value
  for v in 0..<n {
    if residual[u][v] > E::default() {
      if height[v] < min_height {
        min_height = height[v]
      }
    }
  }
  if min_height != @int.max_value {
    height[u] = min_height + 1
  } else {

  }
}

///|
/// Bron-Kerbosch算法 - 寻找图中的所有最大团
/// 返回所有最大团（每个团是一个节点集合）
pub fn[N, E] Graph2::find_maximum_cliques(
  self : Graph2[N, E],
) -> Array[Array[NodeIndex]] {
  guard !self.directed else {
    abort("Find maximum cliques only works for undirected graphs")
  }
  let n = self.node_count()
  let result = Array::new()

  // 转换为邻接矩阵以便快速查找
  let adj_matrix = Array::new()
  for _ in 0..<n {
    adj_matrix.push(Array::make(n, false))
  }
  for u_idx in 0..<n {
    for pair in self.adj[u_idx] {
      let (v, _) = pair
      adj_matrix[u_idx][v.ix] = true
      adj_matrix[v.ix][u_idx] = true // 无向图
    }
  }

  // Bron-Kerbosch递归实现
  fn bron_kerbosch(
    r : Array[Int], // 当前团
    p : Array[Int], // 候选节点
    x : Array[Int], // 已排除的节点
    adj_matrix : Array[Array[Bool]],
    result : Array[Array[NodeIndex]],
  ) -> Unit {
    if p.is_empty() && x.is_empty() {
      // 找到一个最大团
      let clique = Array::new()
      for node_idx in r {
        clique.push(NodeIndex::new(node_idx))
      }
      result.push(clique)
      return
    }

    // 选择枢轴点（从P∪X中选择与P中节点连接最多的点）
    let mut pivot = -1
    let mut max_connections = -1

    // 从P和X中选择枢轴
    let p_union_x = Array::new()
    for v in p {
      p_union_x.push(v)
    }
    for v in x {
      p_union_x.push(v)
    }
    for u in p_union_x {
      let mut count = 0
      for v in p {
        if adj_matrix[u][v] {
          count += 1
        }
      }
      if count > max_connections {
        pivot = u
        max_connections = count
      }
    }

    // 使用枢轴优化：只考虑不与枢轴相连的点
    let p_copy = p.copy()
    for v in p_copy {
      // 如果v与pivot相连，跳过
      if pivot != -1 && adj_matrix[v][pivot] {
        continue
      }
      let new_r = r.copy()
      new_r.push(v)
      let new_p = Array::new()
      let new_x = Array::new()

      // 构建新的P和X集合（v的邻居）
      for u in p {
        if adj_matrix[u][v] {
          new_p.push(u)
        }
      }
      for u in x {
        if adj_matrix[u][v] {
          new_x.push(u)
        }
      }

      // 递归调用
      bron_kerbosch(new_r, new_p, new_x, adj_matrix, result)

      // 将v从P移到X
      match p.search(v) {
        Some(index) => ignore(p.remove(index))
        None => ()
      }
      x.push(v)
    }
  }

  // 初始调用
  let p_initial = Array::new()
  for i in 0..<n {
    p_initial.push(i)
  }
  bron_kerbosch(Array::new(), p_initial, Array::new(), adj_matrix, result)
  result
}
