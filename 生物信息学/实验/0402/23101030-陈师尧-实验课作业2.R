# 定义二叉树结构体函数
createNode <- function(value){
  node <- list(
    value = value,
    left = NULL,
    right = NULL
  )
  return(node)
}

# 插入一个值到二叉树中
insertNode <- function(node, value){
  if(is.null(node)){
    return(createNode(value))
  }
  else if(value < node$value){
    node$left <- insertNode(node$left, value)
  }
  else{
    node$right <- insertNode(node$right, value)
  }
  return(node)
}

# 查找节点
searchNode <- function(node, value){
  if(is.null(node)){
    return(NULL)
  }
  else if(node$value == value){
    return(node)
  }
  else if(node$value > value){
    return(searchNode(node$left, value))
  }
  else{
    return(searchNode(node$right, value))}
}

# 找二叉树中的最小值节点
findMinNode <- function(node){
  if(is.null(node)){
    return(NULL)
  }
  current <- node
  while (!is.null(current$left)) {
    current <- current$left
  }
  return(current)
}

# 删除节点
deleteNode <- function(node, value){
  if(is.null(node)){
    return(NULL)
  }
  else if(value < node$value){
    node$left <- (deleteNode(node$left, value))
    return(node)
  }
  else if(value > node$value){
    node$right <- (deleteNode(node$right, value))
    return(node)
  }
  else{
    if(is.null(node$left)){
      return(node$right)
    }
    else if(is.null(node$right)){
      return(node$left)
    }
    else{
      successor <- findMinNode(node$right)
      node$value <- successor$value
      node$right <- deleteNode(node$right, successor$value)
      return(node)
    }
  }
}

# 二叉树遍历
inorderTraversal <- function(node){
  if(!is.null(node)){
    inorderTraversal(node$left)
    print(node$value)
    inorderTraversal(node$right)
  }
}

preorderTraversal <- function(node){
  if(!is.null(node)){
    print(node$value)
    preorderTraversal(node$left)
    preorderTraversal(node$right)
  }
}

postorderTraversal <- function(node){
  if(!is.null(node)){
    postorderTraversal(node$left)
    postorderTraversal(node$right)
    print(node$value)
  }
}


node <- createNode(50)
node <- insertNode(node, 30)
node <- insertNode(node, 20)
node <- insertNode(node, 40)
node <- insertNode(node, 70)
node <- insertNode(node, 60)
node <- insertNode(node, 80)
cat("原始树的中序遍历:\n")
inorderTraversal(node)
node <- deleteNode(node, 20)
cat("\n删除 20 后的中序遍历:\n")
inorderTraversal(node)
node <- deleteNode(node, 30)
cat("\n删除 30 后的中序遍历:\n")
inorderTraversal(node)
cat("\n删除 50 后的中序遍历:\n")
node <- deleteNode(node, 50)
inorderTraversal(node)