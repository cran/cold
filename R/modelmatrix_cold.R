setMethod("model.mat",
    signature(object = "cold"),
    function (object) 
    {return(object@model.matrix)})