function jointEntropy = CalJointEmtropy(nodeList, idtEntropy, corrEntropy)
    numNodes = length(nodeList);
    m_corrMatrix = inf*ones(numNodes,numNodes);
    % Create partial correlation matrix
    for j=1:numNodes
        node_j = nodeList(j);
        for k=1:numNodes
            node_k = nodeList(k);
            if node_j == node_k
                m_corrMatrix(j,k) = idtEntropy(node_j);
            else
                m_corrMatrix(j,k) = corrEntropy(node_j,node_k);
            end
        end
    end
    % Start calculate joint entropy (greedy method)
    jointEntropy = 0;
    m_candidateList = [];
    for c=1:numNodes
        m_candidateList = [m_candidateList m_corrMatrix(c,c)];
    end
    [m_idtEntropy m_index] = sort(m_candidateList, 'ascend');
    jointEntropy = jointEntropy + m_idtEntropy(1);
    m_corrMatrix(m_index(1),:) = inf;
    for p=2:numNodes
        m_candidateList = m_corrMatrix(:,m_index(1))';
        [m_condEntropy m_index] = sort(m_candidateList, 'ascend');
        jointEntropy = jointEntropy + m_condEntropy(1);
        m_corrMatrix(m_index(1),:) = inf;
    end
end