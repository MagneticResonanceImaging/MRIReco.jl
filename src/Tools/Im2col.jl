#
# rearrange distinct image blocks into columns of a matrix
#
function im2col_distinct(A::Array{T}, blocksize::NTuple{2,Int64}) where T

  nrows = blocksize[1]
  ncols = blocksize[2]
  nelem = nrows*ncols

  # padding for A such that patches can be formed
  row_ext = mod(size(A,1),nrows)
  col_ext = mod(size(A,2),nrows)
  pad_row = (row_ext != 0)*(nrows-row_ext)
  pad_col = (col_ext != 0)*(ncols-col_ext)

  # rearrange matrix
  A1 = zeros(T, size(A,1)+pad_row, size(A,2)+pad_col)
  A1[1:size(A,1),1:size(A,2)] = A

  t1 = reshape( A1,nrows, floor(Int,size(A1,1)/nrows), size(A,2) )
  t2 = reshape( permutedims(t1,[1 3 2]), size(t1,1)*size(t1,3), size(t1,2) )
  t3 = permutedims( reshape( t2, nelem, floor(Int,size(t2,1)/nelem), size(t2,2) ),[1 3 2] )
  res = reshape(t3,nelem,size(t3,2)*size(t3,3))

  return res
end

#
# rearrange overlapping blocks into columns of a matrix
#
function im2col_sliding(A::Array{T}, blocksize::NTuple{2,Int64}) where T
  nrows = blocksize[1]
  ncols = blocksize[2]
  m,n = size(A)

  # start index corresponding to the upper left corner of each block
  start_idx = vec(collect(1:m-nrows+1).+transpose(collect(0:m:m*(n-ncols))))
  # index of a block relative to a starting index
  block_idx = vec(collect(0:nrows-1).+transpose(collect(0:m:m*(ncols-1))))

  res = A[block_idx.+start_idx']
  return res
end

#
# rearrange columns of a matrix into blocks of an image
#--------------------------------------------------------
# size(A) should not be larger then (blocksize[1]*blocksize[2], matsize[1]*matsize[2]).
# otherwise the bottom (right) lines (columns) will be cut.
# matsize should be divisble by blocksize.
#
function col2im_distinct(A::Array{T}, blocksize::NTuple{2,Int64}, matsize::NTuple{2,Int64}) where T

  if mod(matsize[1],blocksize[1]) != 0 || mod(matsize[2],blocksize[2]) != 0
    error("matsize should be divisible by blocksize")
  end

  blockrows = blocksize[1]
  blockcols = blocksize[2]
  matrows = matsize[1]
  matcols = matsize[2]
  mblock = floor(Int,matrows/blockrows) # number of blocks per row
  nblock = floor(Int,matcols/blockcols) # number of blocks per column

  # padding for A such that patches can be formed and arranged into a matrix of
  # adequate size
  row_ext = blockrows*blockcols-size(A,1)
  col_ext = mblock*nblock-size(A,2)
  pad_row = (row_ext > 0 )*row_ext
  pad_col = (col_ext > 0 )*col_ext

  A1 = zeros(T, size(A,1)+pad_row, size(A,2)+pad_col)
  A1[1:blockrows*blockcols, 1:mblock*nblock] = A[1:blockrows*blockcols, 1:mblock*nblock]

  # rearrange matrix
  t1 = reshape( A1, blockrows,blockcols,mblock*nblock )
  t2 = reshape( permutedims(t1,[1 3 2]), matrows,nblock,blockcols )
  res = reshape( permutedims(t2,[1 3 2]), matrows,matcols)

end

#
# rearrange columns of a matrix into overlapping blocks of an image
#
function col2im_sliding(A::Array{T}, blocksize::NTuple{2,Int64}, matsize::NTuple{2,Int64}) where T
  nrows, ncols = blocksize
  m,n = matsize
  blocksperrow = m-nrows+1
  blockspercol = n-ncols+1

  res = zeros(T,matsize)
  res[1:nrows-1,1:blockspercol] = A[1:nrows-1,1:blocksperrow:end]
  res[nrows:m,1:blockspercol] = reshape(A[nrows,:],m-nrows+1,blockspercol)
  res[1:nrows, blockspercol+1:end] = reshape(A[nrows+1:end,(blockspercol-1)*blocksperrow+1], nrows, ncols-1)
  res[nrows+1:m, blockspercol+1:end] = copy(transpose(A[2*nrows:nrows:end,(blockspercol-1)*blocksperrow+2:end]))

  return res
end

function test_im2col()
  a = reshape(collect(1:24),4,6)
  blocksize = (2,3)

  b = im2col_distinct(a, blocksize)
  return a,b
end

function test_col2im()
  a = reshape(collect(1:24),4,6)
  blocksize = (2,3)

  b = im2col_distinct(a, blocksize)
  c = col2im_distinct(b,blocksize,(4,6))
  return a,b,c
end
