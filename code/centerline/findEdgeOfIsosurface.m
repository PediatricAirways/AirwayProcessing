function vertexOnEdge = findEdgeOfIsosurface( vertex, szf, face )

edges = zeros( szf*3, 2 );
for iI = 1:szf
    edges( ( iI - 1 ) * 3 + 1, 1 ) = min( face( iI, 1 ), face( iI, 2 ) );
    edges( ( iI - 1 ) * 3 + 1, 2 ) = max( face( iI, 1 ), face( iI, 2 ) );
    edges( ( iI - 1 ) * 3 + 2, 1 ) = min( face( iI, 2 ), face( iI, 3 ) );
    edges( ( iI - 1 ) * 3 + 2, 2 ) = max( face( iI, 2 ), face( iI, 3 ) );
    edges( ( iI - 1 ) * 3 + 3, 1 ) = min( face( iI, 3 ), face( iI, 1 ) );
    edges( ( iI - 1 ) * 3 + 3, 2 ) = max( face( iI, 3 ), face( iI, 1 ) );
end

edgesSorted = sortrows( edges, 2 );
edgesSorted = sortrows( edgesSorted, 1 );

count = 0;
iBegin = 1;
while iBegin <= length( edgesSorted )
    iEnd = iBegin+1;
    while iEnd <= length( edgesSorted )
        if edgesSorted( iEnd, 1 ) == edgesSorted( iBegin, 1 ) & edgesSorted( iEnd, 2 ) == edgesSorted( iBegin, 2 )
            iEnd = iEnd + 1;
        else
            break;
        end
    end
    if iEnd - iBegin == 1
        count = count + 1;
        vertexIdOnEdge( count, : ) = edgesSorted( iBegin, : ); 
    end
    iBegin = iEnd;
end

vertexIdOnEdge = sort( reshape( vertexIdOnEdge, count*2, 1 ) );
count = 0;
iBegin = 1;
while iBegin <= length( vertexIdOnEdge )
    count = count + 1;
    vertexOnEdge( count, : ) = vertex( vertexIdOnEdge( iBegin ), : );
    iEnd = iBegin + 1;
    while iEnd <= length( vertexIdOnEdge )
        if vertexIdOnEdge( iEnd ) == vertexIdOnEdge( iBegin )
            iEnd = iEnd + 1;
        else
            break;
        end
    end
    iBegin = iEnd;
end
