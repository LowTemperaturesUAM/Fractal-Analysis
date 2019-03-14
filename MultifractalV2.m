%Function MultifractalV2 with:
%Inputs -     Matriz: Square and 2 exponent size Matrix
%             Q:      Array of the exponents of Q, tipically Q is symetric
%                     Qmax = 10 and 1000 points
%             and 

%Outputs -    A:      Matrix of sum(mu*P) with dependance with eps and Q
%             eps:    Array of the sizes of the box, also exponents of 2 
%             TAU:    Matrix of the different values of tau with a
%                     dependance in eps and Q
%             IQ:     Matrix of sum(P^q) with dependance in Q and eps (I
%                     dont use it for the rest


function [A, epsilon, TAU, IQ] = MultifractalV2( Matriz, Q )


                    % Pad the image with background pixels so that its dimensions are a power of 2.
     maxDim = max(size(Matriz));
                    %     newDimSize = 2^ceil(log2(maxDim));
                    %     rowPad = newDimSize - size(I, 1);
                    %     colPad = newDimSize - size(I, 2);
                    %     I = padarray(I, [rowPad, colPad], 'post');

    %Predefino los dos vectores salidas con cero
    %Predefine the output arrays with ceros
    boxCounts = zeros(1, ceil(log2(maxDim)));
    resolutions = zeros(1, ceil(log2(maxDim)));
    
    boxSize = size(Matriz, 1);
    boxesPerDim = 1;
    idx = 0;
    
    %BUcle donde el tamaño de la caja se va haciendo mas pequea
    %Looping and decreasing the box size
    while boxSize >= 1
        boxCount = 0;
        numberperBox = 0;
        idx = idx + 1;
        %Bucles donde van moviendo la caja. Primero por filas, y luego
        %cambia a la siguiente columna y continua
        %This loop moves the box. First in rows and then in columns.
        for boxRow = 1:boxesPerDim
            for boxCol = 1:boxesPerDim
                
             
                minRow = (boxRow - 1) * boxSize + 1;
                maxRow = boxRow * boxSize;
                minCol = (boxCol - 1) * boxSize + 1;
                maxCol = boxCol * boxSize;
                
            %Comprobacion del movimiento de las cajas 
            %Uncoment to check how the box is moving.
%                 figure(11);
%                 imshow(I)
%                 hold on
%                 rectangle('Position', [minRow, minCol, maxRow-minRow, maxCol-minCol],'EdgeColor', 'w')
                
                
                %Conteo del numero de cajas
                %Number of boxes at same size counter
                boxCount = boxCount + 1;
                %Vedctor que te guarda el numero de puntos encontrados en
                %cada caja, por tanto su dmension es el numero de cajas
                
                %Array for saving the number of found pixel in each box at the same size
                numberperBox(boxCount) = sum(sum(Matriz(minCol:maxCol,minRow:maxRow)));
                
                
          
          
            end

        end
        %Calculo del numero de puntos por caja entre el numero de puntos
        %totales.
        %Calculate the probability of finding the points in the box.
        ProbabilityperBox = numberperBox./sum(numberperBox);
        ProbabilityperBox = nonzeros(ProbabilityperBox);
        boxSize
        %Numero de cajas que contienen mas de 1 punto
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
        
            
             
             
        
        %Disminuyo por la mitad el tamaño de cajas
        %Decrease by half the box size.
        boxesPerDim = boxesPerDim * 2;
        boxSize = boxSize / 2;
    end


   
     
end
