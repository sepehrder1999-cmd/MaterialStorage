function i =RouletteWheelSelection(P)

    r=rand;

    c=cumsum(P);

    i=find(r<=c,1,'first');
    %shomare andisi ro mide ke baraye avalin bar in shart dorost mishe
end