function clattice=clattice(resolution,numx,numy,radius,spacing,offset) % Generates a lattice of circles to be used as a scene.
n=(1:resolution)∗2∗pi/resolution; % resolution= # of segments forming each circle.
x = radius∗cos(n)+offset(1);  % offset(1,2)= (x,y) coordinates.
y = radius∗sin(n)+offset(2);  % numx, numy = lattice size.
clattice=[];  
for i = 1:numx
  for j = 1:numy
    if mod(j,2)==0
      clattice=[clattice;formatobject([x"+spacing∗(i−1)+spacing/2,y"+spacing∗(j−1)])];
    else
      clattice=[clattice;formatobject([x"+spacing∗(i−1),y"+spacing∗(j−1)])];
    end
  end
end
