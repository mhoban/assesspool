// generate combinations of list items
def combn(list,m) {
    def n = list.size()
    m == 0 ?
        [[]] :
        (0..(n-m)).inject([]) { newlist, k ->
            def sublist = (k+1 == n) ? [] : list[(k+1)..<n]
            newlist += combn(sublist,m-1).collect { [list[k]] + it }
        }
}
