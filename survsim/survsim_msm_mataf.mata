mata:
function hf${i}(tnodes,expxb,tdexb,hvars,lt,time0) 
{
        real matrix hf
	hf = ${matahazard${i}}
	return(hf)
}
end
