% opl2phase.m computes the phase at the end of a beam's path.
function phase=opl2phase(OPL,lamda) 
phase=mod(OPL,lamda)∗2∗pi/lamda; 
end
