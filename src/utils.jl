#returns date and time
dtnow() = Dates.format(now(), "Yud_HHhMM\\mSS\\s")
function println_timestamp(io::IO, msg)
    println(string(dtnow(),": ", msg))
    flush(io)
end
println_timestamp(msg) = println_timestamp(stdout, msg)


#sorting functions for rank_spus
function InsertionSort!(A::AbstractArray{T, 1}, order::AbstractArray{Int64, 1}, ii=1, jj=length(A)) where {T<:Real}
    for i = ii+1 : jj
        j = i - 1
        temp  = A[i]
        itemp = order[i]
        while true
            if j == ii-1
                break
            end
            if A[j] <= temp
                break
            end
            A[j+1] = A[j]
            order[j+1] = order[j]
            j -= 1
        end
        A[j+1] = temp
        order[j+1] = itemp
    end  # i
    return A
end # function InsertionSort!

function quicksort!(A::AbstractArray{T, 1}, order::AbstractArray{Int64, 1}, i=1, j=length(A)) where {T<:Real}
    if j > i
      if  j - i <= 10
        # Insertion sort for small groups is faster than Quicksort
        InsertionSort!(A,order, i,j)
        return A
      end
      pivot = A[ div(i+j,2) ]
      left, right = i, j
      while left <= right
        while A[left] < pivot
            left += 1
        end
        while A[right] > pivot
            right -= 1
        end
        if left <= right
            A[left], A[right] = A[right], A[left]
            order[left], order[right] = order[right], order[left]
            left += 1
            right -= 1
        end
      end  # left <= right
      quicksort!(A,order, i, right)
      quicksort!(A,order, left, j)
    end  # j > i
    return A
end # function quicksort!
