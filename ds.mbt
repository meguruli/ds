///|
pub struct NodeIndex {
  ix : Int
} derive(Eq, Compare, Show, Hash)

pub fn NodeIndex::new(ix : Int) -> NodeIndex {
  NodeIndex::{ ix, }
}

///| 基于邻接表的简单图实现（用于作业扩展）
pub struct Graph2[N, E] {
  directed : Bool
  adj : Array[Array[(NodeIndex, E)]]
  node_weights : Array[N]
  edges : Array[(NodeIndex, NodeIndex, E)]
} derive(Show)

pub fn[N : Default, E : Default] Graph2::new(directed~ : Bool = true) -> Graph2[N, E] {
  Graph2::{
    directed,
    adj: Array::new(),
    node_weights: Array::new(),
    edges: Array::new(),
  }
}

pub fn[N, E : Default] Graph2::node_count(self : Graph2[N, E]) -> Int {
  self.adj.length()
}

pub fn[N : Default, E : Default] Graph2::add_node(self : Graph2[N, E], weight? : N) -> NodeIndex {
  let w = weight.unwrap_or(N::default())
  let idx = self.adj.length()
  self.node_weights.push(w)
  self.adj.push(Array::new())
  NodeIndex::new(idx)
}

pub fn[N, E : Default] Graph2::add_edge(self : Graph2[N, E], a : NodeIndex, b : NodeIndex, weight? : E) -> Unit {
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

pub fn[N, E] Graph2::neighbors(self : Graph2[N, E], a : NodeIndex) -> Array[NodeIndex] {
  let res = Array::new()
  for pair in self.adj[a.ix] {
    res.push(pair.0)
  }
  res
}

///| 无向连通分量（返回每个连通分量的节点列表）
pub fn[N, E] Graph2::connected_components(self : Graph2[N, E]) -> Array[Array[NodeIndex]] {
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

///| Kosaraju 强连通分量（适用于有向图）
pub fn[N, E] Graph2::strongly_connected_components(self : Graph2[N, E]) -> Array[Array[NodeIndex]] {
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
        if vis2[u.ix] { continue }
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
    if jj == 0 { break }
  }
  sccs
}

///| Dijkstra（仅针对 Int 权重实现）
pub fn Graph2::dijkstra_int(self : Graph2[Int, Int], a : NodeIndex) -> Array[Int] {
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
      if vis[v.ix] { continue }
      let nd = dist[u] + w
      if nd < dist[v.ix] {
        dist[v.ix] = nd
      }
    }
  }
  dist
}

///| Kruskal 最小生成树（仅针对 Int 权重的无向图）
pub fn[N] Graph2::kruskal_mst(self : Graph2[N, Int]) -> Array[(NodeIndex, NodeIndex, Int)] {
  let res = Array::new()
  let n = self.adj.length()
  // collect edges (undirected) ensuring each undirected edge counted once
  let es = Array::new()
  for e in self.edges {
    let (u, v, w) = e
    if self.directed {
      es.push((u, v, w))
    } else {
      if u.ix <= v.ix {
        es.push((u, v, w))
      } else {
        es.push((v, u, w))
      }
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
    if ru == rv { continue }
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

///| Traversal helpers adapted from NetworkX (DFS, BFS, Topo)
struct Dfs {
  stack : Array[NodeIndex]
  discovered : VisitMap
}

pub fn[N, E : Default] Dfs::new(graph : Graph2[N, E], start : NodeIndex) -> Dfs {
  Dfs::{ stack: [start], discovered: graph.visit_map() }
}

struct VisitMap {
  visit : Array[Bool]
}

pub fn VisitMap::new() -> VisitMap {
  VisitMap::{ visit: Array::new() }
}

fn VisitMap::contain(self : VisitMap, nx : Int) -> Bool {
  self.visit[nx]
}

impl VM for VisitMap with visit(self, node : NodeIndex) -> Bool {
  let nx = node.ix
  let res = self.contain(nx)
  self.visit[nx] = true
  not(res)
}

impl VM for VisitMap with is_visited(self, node : NodeIndex) -> Bool {
  self.contain(node.ix)
}

priv trait VM {
  visit(Self, NodeIndex) -> Bool
  is_visited(Self, NodeIndex) -> Bool
}

priv trait Visitable {
  visit_map(Self) -> VisitMap
}

impl[N, E : Default] Visitable for Graph2[N, E] with visit_map(self) -> VisitMap {
  VisitMap::{ visit: Array::make(self.node_count(), false) }
}

fn VisitMap::clear(self : VisitMap) -> Unit {
  self.visit.clear()
}

pub fn Dfs::from_parts(stack : Array[NodeIndex], discovered : VisitMap) -> Dfs {
  Dfs::{ stack, discovered }
}

pub fn Dfs::reset(self : Dfs) -> Unit {
  self.stack.clear()
  self.discovered.clear()
}

pub fn Dfs::move_to(self : Dfs, node : NodeIndex) -> Unit {
  self.stack.clear()
  self.stack.push(node)
}

pub fn[N, E : Default] Dfs::next(self : Dfs, graph : Graph2[N, E]) -> NodeIndex? {
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

struct DfsPostOrder {
  stack : Array[NodeIndex]
  discovered : VisitMap
  finished : VisitMap
}

pub fn[N, E : Default] DfsPostOrder::new(
  graph : Graph2[N, E],
  start : NodeIndex
) -> DfsPostOrder {
  DfsPostOrder::{
    stack: [start],
    discovered: graph.visit_map(),
    finished: graph.visit_map(),
  }
}

pub fn DfsPostOrder::reset(self : DfsPostOrder) -> Unit {
  self.stack.clear()
  self.discovered.clear()
  self.finished.clear()
}

pub fn DfsPostOrder::move_to(self : DfsPostOrder, node : NodeIndex) -> Unit {
  self.stack.clear()
  self.stack.push(node)
}

pub fn[N, E : Default] DfsPostOrder::next(
  self : DfsPostOrder,
  graph : Graph2[N, E]
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

struct Bfs {
  stack : @deque.Deque[NodeIndex]
  discovered : VisitMap
}

pub fn[N, E : Default] Bfs::new(graph : Graph2[N, E], start : NodeIndex) -> Bfs {
  let stack = @deque.new()
  stack.push_front(start)
  let discovered = graph.visit_map()
  discovered.visit(start) |> ignore
  Bfs::{ stack, discovered }
}

pub fn[N, E : Default] Bfs::next(self : Bfs, graph : Graph2[N, E]) -> NodeIndex? {
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

struct Topo {
  tovisit : Array[NodeIndex]
  ordered : VisitMap
  degrees : Array[Int]
}

pub fn[N, E : Default] Topo::new(graph : Graph2[N, E]) -> Topo {
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

pub fn[N, E : Default] Topo::empty(graph : Graph2[N, E]) -> Topo {
  Topo::{
    tovisit: Array::new(),
    ordered: graph.visit_map(),
    degrees: Array::new(),
  }
}

pub fn[N, E : Default] Topo::next(self : Topo, graph : Graph2[N, E]) -> NodeIndex? {
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

fn[N, E : Default] get_indegree(self : Graph2[N, E]) -> Array[Int] {
  let deg = Array::make(self.node_count(), 0)
  for i in 0..<self.node_count() {
    for pair in self.adj[i] {
      let (v, _) = pair
      deg[v.ix] += 1
    }
  }
  deg
}

///| Generic Dijkstra adapted for Graph2
pub fn[N, E : Compare + Add + Default] Graph2::dijkstra(
  self : Graph2[N, E],
  a : NodeIndex
) -> @hashmap.HashMap[NodeIndex, E] {
  let n = self.adj.length()
  guard a.ix >= 0 && a.ix < n else {
    abort("Node not found")
  }
  let vis = Array::make(n, false)
  let dist = @hashmap.new()
  let q : @priority_queue.T[(E, NodeIndex)] = @priority_queue.new()
  q.push((E::default(), a))
  while not(q.is_empty()) {
    let (d, u) = q.pop().unwrap()
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
        q.push((new_dist, v))
      }
    }
  }
  dist
}

///| LCA adapted for Graph2 (removed cycle check to match Graph2 capabilities)
struct Lca[N, E] {
  graph : Graph2[N, E]
  deep : Array[Int]
  father : Array[Array[Int]]
  mut root : NodeIndex
}

fn[N, E : Default] Lca::update_lca_info_pre(
  self : Lca[N, E],
  graph : Graph2[N, E],
  u : NodeIndex,
  f : NodeIndex,
  d : Int
) -> Unit {
  self.deep[u.ix] = d
  for v in graph.neighbors(u) {
    if v != f {
      self.father[v.ix][0] = u.ix
      self.update_lca_info_pre(graph, v, u, d + 1)
    }
  }
}

fn[N, E : Default] Lca::update_lca_info(
  self : Lca[N, E],
  graph : Graph2[N, E],
  u : NodeIndex,
  f : NodeIndex,
  d : Int
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

pub fn[N, E : Default] Lca::new(
  graph : Graph2[N, E],
  root : NodeIndex
) -> Lca[N, E] {
  let n = graph.node_count()
  let deep = Array::make(n + 1, 0)
  let top = @math.log2((n + 1).to_double()).to_int()
  let father = Array::makei(n + 1, fn(_) { Array::make(top + 1, 0) })
  let lca = Lca::{ graph, deep, father, root }
  lca.update_lca_info(lca.graph, lca.root, NodeIndex::new(-1), 1)
  lca
}

pub fn[N, E : Default] Lca::query(
  self : Lca[N, E],
  u : NodeIndex,
  v : NodeIndex
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

pub fn[N, E : Default] Lca::reset(self : Lca[N, E], root? : NodeIndex) -> Unit {
  if root is Some(r) {
    self.root = r
  }
  self.update_lca_info(self.graph, self.root, NodeIndex::new(-1), 1)
}

///| Utility tuple helpers copied from NetworkX
fn[T] get_tuple(tuple : (T, T), index : Int) -> T {
  if index == 0 {
    return tuple.0
  }
  tuple.1
}

///|
fn[T] set_tuple(tuple : (T, T), index : Int, value : T) -> (T, T) {
  if index == 0 {
    return (value, tuple.1)
  }
  (tuple.0, value)
}

///| Error types used by network utilities
priv suberror NodeIxLimit String derive(Show)

///|
priv suberror NodeOutBounds String derive(Show)

///|
priv suberror NodeNotFound String derive(Show)

///|
suberror CycleDetected String derive(Show)

pub fn[N, E] Graph2::find_edge(self : Graph2[N, E], a : NodeIndex, b : NodeIndex) -> Int? {
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

pub fn[N, E : Eq] Graph2::remove_edge_nodes(self : Graph2[N, E], a : NodeIndex, b : NodeIndex) -> E? {
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

pub fn[N, E : Eq + Default] Graph2::remove_node(self : Graph2[N, E], a : NodeIndex) -> N? {
  let n = self.node_count()
  guard a.ix >= 0 && a.ix < n else { return None }

  // remove all incident edges (works for directed and undirected)
  while true {
    let mut found = false
    let mut ei = 0
    while ei < self.edges.length() {
      let (u, v, _) = self.edges[ei]
      if u == a || v == a {
        // remove this edge (will update adj and edges)
        self.remove_edge_nodes(u, v) |> ignore
        found = true
        break
      }
      ei = ei + 1
    }
    if not(found) { break }
  }

  // grab node weight to return
  let weight = self.node_weights[a.ix]

  // swap-remove the node
  let last_index = self.node_count() - 1
  if a.ix != last_index {
    // move last node into position a
    self.node_weights.swap(a.ix, last_index)
    self.adj.swap(a.ix, last_index)

    let old_idx = NodeIndex::new(last_index)
    let new_idx = NodeIndex::new(a.ix)

    // update adjacency lists: replace any reference to old_idx with new_idx
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

    // update edges array: replace node occurrences
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

  // remove last slot
  ignore(self.node_weights.remove(last_index))
  ignore(self.adj.remove(last_index))
  Some(weight)
}

pub fn[N, E] Graph2::edge_count(self : Graph2[N, E]) -> Int {
  self.edges.length()
}

pub fn[N, E] Graph2::find_edge_undirected(self : Graph2[N, E], a : NodeIndex, b : NodeIndex) -> Int? {
  // Graph2::find_edge already handles undirected graphs, so reuse it
  self.find_edge(a, b)
}

pub fn[N, E] Graph2::change_edge_links(
  self : Graph2[N, E],
  edge_node : (NodeIndex, NodeIndex),
  e : Int,
  edge_next : (Int, Int)
) -> Unit {
  // Graph2 uses explicit adjacency lists (self.adj) and does not maintain
  // per-edge successor links like `Graph`, so there is nothing to update here.
}

pub fn[N, E] Graph2::remove_edge_adjust_indices(self : Graph2[N, E], e : Int) -> E? {
  let last_index = self.edges.length() - 1
  if e != last_index {
    self.edges.swap(e, last_index)
  }
  let (u, v, w) = self.edges.remove(last_index)
  Some(w)
}

