corner <- function(input.mat){
    if(nrow(input.mat) > 4){
        row.out = 5
    } else { row.out = nrow(input.mat)}

    if(ncol(input.mat) > 7){
        col.out = 8
    } else { col.out = ncol(input.mat)}

    print(input.mat[1:row.out, 1:col.out])
}
