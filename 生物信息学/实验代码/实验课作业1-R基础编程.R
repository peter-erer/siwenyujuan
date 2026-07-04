# 1.定义一个包含从1到100的向量x
x <- 1:100

# 2.将向量x中的偶数定义成向量y
y <- x [x %% 2 ==0]

# 3.将向量x中的奇数定义成向量z
z <- x [x %% 2 ==1]

# 4.创建一个50x2的矩阵m，该矩阵包含y，z，并将列分别命名成Even和Odd
m <- cbind(Even <- y, Odd <- z)

# 5.利用x，创建一个10x10的矩阵n，并将列分别命名成C1，C2，…，C10，将行命名成R1，R2，…，R10
n <- matrix(x, 10, 10)
colnames(n) <- paste0("C",1:10)
rownames(n) <- paste0("R",1:10)
n

# 6.创建一个名为list_data的list，该list中包含x，y，z，m，n，并分别命名成vector_x，vector_y, vector_z，matrix_m，matrix_y，让后输出list中的matrix_y
list_data <- list(vector_x = x, 
                  vector_y = y, 
                  vector_z = z, 
                  matrix_m = m, 
                  matrix_y = n)
list_data$matrix_y

# 7. 将x中的数分成5给等级A，B，C，D，E（A>= 90；80<=B<90; 70<=C<80；60<=D<70；E<60）
grade_fun <- function(vec){
  res <- character(length(vec))
  res[vec >= 90] <- "A"
  res[vec >= 80 & vec < 90] <- "B"
  res[vec >= 70 & vec < 80] <- "C"
  res[vec >= 60 & vec < 70] <- "D"
  res[vec < 60] <- "E"
  return(res)
}
grade_fun(x)

# 8.统计y中A，B，C，D，E的数目
count_fun <- function(vec){
  vec_factor <- factor(grade_fun(vec),levels = c("A","B","C","D","E"))
  res <- table(vec_factor)
  return(res)
}
count_fun(y)

# 9.统计z中A，B，C，D，E的数目
count_fun(z)

# 10.编写函数，利用向量y和z分别按公式计算
calc_E <- function(v){
  E <- 1
  for (i in 1:length(v)) {
    E <- E + v[i]^i/factorial(i)
  }
  return(E)
}
Ey <- calc_E(y)
Ez <- calc_E(z)
dif <- Ey - Ez
if (dif > 0) {
  print("Ey > Ez")
}else if (dif < 0) {
  print("Ey < Ez")
}else{
  print("Ey = Ez")
}
abs_dif <- abs(dif)
abs_dif