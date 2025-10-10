///|
pub struct NodeIndex {
  ix : Int
} derive(Eq, Compare, Show, Hash)

///|
pub fn NodeIndex::new(ix : Int) -> NodeIndex {
  NodeIndex::{ ix, }
}

///|
/// 基于邻接表的简单图实现
pub struct Graph2[N, E] {
  directed : Bool
  adj : Array[Array[(NodeIndex, E)]]
  node_weights : Array[N]
  edges : Array[(NodeIndex, NodeIndex, E)]
} derive(Show)

///|
pub fn[N, E] Graph2::new(directed? : Bool = true) -> Graph2[N, E] {
  Graph2::{
    directed,
    adj: Array::new(),
    node_weights: Array::new(),
    edges: Array::new(),
  }
}
