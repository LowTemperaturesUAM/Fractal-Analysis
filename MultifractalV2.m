%Function MultifractalV2 with:
%Inputs -     Matriz: Matrix has to be square and with a number of lines being an exponent of 2
%             Q:      Array of exponents of Q. Tipically Q is symmetric with respect to zero
%                     Qmax = 10 and 1000 points
%             and 

%Outputs -    A:      Matrix of sum(mu*P) which depends on eps and Q
%             eps:    Array of the sizes of the box. This will always be exponents of 2 
%             TAU:    Matrix of the different values of tau, depends on eps and Q
%             IQ:     Matrix of sum(P^q), which depends on eps and Q (I dont use it for the rest)


function [A, epsilon, TAU, IQ] = MultifractalV2( Matriz, Q )


                    % Pad the image with background pixels so that its dimensions are a power of 2.
     maxDim = max(size(Matriz));
                    %     newDimSize = 2^ceil(log2(maxDim));
                    %     rowPad = newDimSize - size(I, 1);
                    %     colPad = newDimSize - size(I, 2);
                    %     I = padarray(I, [rowPad, colPad], 'post');

    %Predefine the output arrays with ceros
    boxCounts = zeros(1, ceil(log2(maxDim)));
    resolutions = zeros(1, ceil(log2(maxDim)));
    
    boxSize = size(Matriz, 1);
    boxesPerDim = 1;
    idx = 0;
    
    %Loop and decreasing the box size
    while boxSize >= 1
        boxCount = 0;
        numberperBox = 0;
        idx = idx + 1;
        %This loop moves the box. First in rows and then in columns.
        for boxRow = 1:boxesPerDim
            for boxCol = 1:boxesPerDim
                
             
                minRow = (boxRow - 1) * boxSize + 1;
                maxRow = boxRow * boxSize;
                minCol = (boxCol - 1) * boxSize + 1;
                maxCol = boxCol * boxSize;
                
            %Uncomment to check how the box is moving.
%                 figure(11);
%                 imshow(I)
%                 hold on
%                 rectangle('Position', [minRow, minCol, maxRow-minRow, maxCol-minCol],'EdgeColor', 'w')
                
                
                %Number of boxes at same size counter
                boxCount = boxCount + 1;
                %Array for saving the number of found pixel in each box at the same size
                numberperBox(boxCount) = sum(sum(Matriz(minCol:maxCol,minRow:maxRow)));
                
                
          
          
            end

        end
        %Calculate the probability of finding the points in the box.
        ProbabilityperBox = numberperBox./sum(numberperBox);
        ProbabilityperBox = nonzeros(ProbabilityperBox);
        boxSize
        %Number of boxes with more than 1 point.
        NumberOfmore0Pixel = sum(numberperBox > 0);
        
        ProbabilityQ = zeros(length(Q), length(ProbabilityperBox));
        for i =1:length(Q)
     %      I(i) = sum(ProbabilityperBox.^(Q(i)))./NumberOfmore0Pixel;
           I(i) = sum(ProbabilityperBox.^(Q(i)));
           
          
           tau1(i) = sum(ProbabilityperBox.^(Q(i)-1))./NumberOfmore0Pixel;

           ProbabilityQ(i,:) = ProbabilityperBox.^Q(i);
        end
  
        MU = ProbabilityQ./I';
        IQ(idx,:) = I;
       % A(idx,:) = sum(MU.*log(ProbabilityperBox'),2);
        A(idx,:) = sum(MU.*ProbabilityperBox',2);

        
        %MU = ProbabilityperBox.^Q
       
        TAU(idx,:) = tau1;
        epsilon(idx) = boxSize;
        
            
             
             
        
        %Decrease by half the box size.
        boxesPerDim = boxesPerDim * 2;
        boxSize = boxSize / 2;
    end


   
     
end
