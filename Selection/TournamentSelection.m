function i=TournamentSelection(pop,m)

    nPop=numel(pop);

    %majmooe oonai ke entekhab shodan
    S=randsample(nPop,m);

    spop=pop(S);

    scosts=[spop.Cost];

    % ma ba khode meghdare min kari nadarim ba mahalesh kar darim ke andise
    % j hast
    [~, j]=min(scosts);

    i=S(j);

end