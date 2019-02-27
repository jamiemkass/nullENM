#' @export
nullSDMResults <- setClass("nullSDMResults",
                          slots=c(model='character',
                                  model.args='list',
                                  eval.type='character',
                                  occs='data.frame',
                                  occs.grp='numeric',
                                  bg='data.frame',
                                  bg.grp='numeric',
                                  no.iter='numeric',
                                  all.stats='data.frame',
                                  null.stats='data.frame',
                                  null.stats.iters='data.frame'))
#' @export
setMethod("show",
          signature="nullSDMResults",
          definition=function(object) {
            cat("An object of class: ", class(object), "\n", sep = "")
            cat(" ",  "Model: ", object@model, '\n',
                " ", "Occurrence/background points: ",
                paste0(nrow(object@occs), '/', nrow(object@bg)), '\n',
                " ",  "Evaluation type: ", object@eval.type, '\n',
                " ",  "Number of iterations: ", object@no.iter, '\n', sep = "")
            invisible(NULL)
          })
