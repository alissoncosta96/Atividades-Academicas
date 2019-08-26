#DUPLA ALISSON COSTA 


############## ARMAZENAMENTO DE DADOS - INICIO ################
C = input('Digite a Função Objetivo na forma vetorial [X1 X2 ... XN]: ');#ENTRADA 'função objetivo'
disp('');
m = input('Digite a quantidade de restrições: '); #ENTRADA 'quantidade restrições'

for i= 1:m
    Ab(i,:)=[input('Digite a restrição na forma vetorial: ')]; #ARMAZENA AS restrições JUNTAMENTE 
                                                               #COM OS 'termos após a igualdade'
end

nAb = size(Ab, 2); #QUANTIDADE DE COLUNAS DE Ab

for i=1:nAb-1   
     A(:,i) = Ab(:,i); #RETIRA A ultima coluna 'b' 
                       #DA MATRIZ Ab SEPARANDO 'A' e 'b' - MATRIZ 'A'
end

b = Ab(:, nAb); #ARMAZENA A MATRIZ 'b' 

n = size(A,2); #QUANTIDADE DE COLUNAS DE 'A'
varIndep = n - m; #QUANTIDADE DE VARIÁVEIS INDEPENDENTES QUE PODEM SER FIXADAS COMO ZERO

k=0;    
for j=1:n 
    k = k+1;
    indiceVar(1, j) = k;#ATRIBUI 'índices' A TODAS AS 'colunas (variáveis)'
end  

Xb = zeros(1,n-varIndep);#INICIALIZANDO AS 'variáveis básicas' COM 'zero'
############# ARMAZENAMENTO DE DADOS - FIM #################
 

#INÍCIO DO SIMPLEX
g = 1;
contAux = 1; 
while g==1
  p = 1;
    while p==1   
        if contAux==1#VERIFICA SE os índices B e N 
            if isempty(find(Xb<0))# VERICA SE TODOS OS ELEMENTOS 
                                  # DE 'Xb' SÃO POSITIVOS (FACTIBILIDADE)
                                  
                indN = indiceVar(1:varIndep);#DEFINE OS 'índices das variáveis básicas'
                                             #PARTIÇÃO NÃO BÁSICA
                                             
                indB = indiceVar(varIndep+1:n);#DEFINE OS 'índices das variáveis não básicas'
                                               #PARTIÇÃO BÁSICA
                
            else #CASO Xb 'não seja factível'
                indN = indiceVar(n-varIndep+1:n); #DEFINE 'novos índices das variáveis básicas'
                indB = indiceVar(1:n-varIndep);   #DEFINE 'novos índices das variáveis não básicas'                
            end  
        end 
        
        B = A(:, indB); #MATRIZ 'B'
        N = A(:, indN); #MATRIZ 'N'
            
        
       Xb=inv(B)*b;   
       if isempty(find(Xb<0)) #VERIFICA A 'factibilidade' DA SOLUÇÃO BÁSICA     
           p=0;
       end    
    endwhile

    Cb = C(indB);#'função objetivo' ESCRITA NOS    
                 #'índices das variáveis básicas'
                 
    Cn = C(indN);#'função objetivo' ESCRITA NOS 
                 #'índices das variáveis não básicas'

    mSimplex =  Cb*inv(B); #'Multiplicador Simplex': NECESSÁRIO PARA
                           # CALCULAR OS 'Custos Reduzidos'
    
    
    q=1;   
    
    while q<=numel(indN)
        auxRed = Cn(1,q) - ((mSimplex)*(N(:, q))); #  CALCULA OS 'Custos Reduzidos'
        cRed(1,q) = auxRed; # ARMAZENA OS 'Custos Reduzidos' NO VETOR
        auxRed = [];
        q++;     
    endwhile
    
    [menorCR, indK] = min(cRed);# ARMAZENA O 'menor custo reduzido' COM SEU 'índice'
                                # DETERMINAÇÃO DO ÍNDICE A ENTRAR NA BASE  
                                
    entraNaBase =  indN(1,indK);# ÍNDICE QUE VAI ENTRAR NA 'base'
    
    y = inv(B)*N(:,indK);#  DIREÇÃO SIMPLEX
    
    if isempty(find(y>0))#  VERIFICA SE TODOS O ELEMENTOS 
                         #  DA 'Direção Simplex' SÃO NEGATIVOS (SOLUÇÃO ILIMITADA)
          disp('-------------------------------------------------------------------------');
          disp("A SOLUÇÃO ÓTIMA É ILIMITADA COM A FUNÇÃO CUSTO AUMENTANDO PARA O INFINITO");
          disp('-------------------------------------------------------------------------');
          break;
    end
    
    if isempty(find(cRed<0))# VERIFICA SE TODOS OS ELEMENTOS
                            # DO 'Custo Reduzido' SÃO POSITIVOS (SOLUÇÃO ÓTIMA) 
        g=0;
        z = Cb*Xb; # CUSTO ÓTIMO
        xSol = zeros(1, n);
        xSol(indB) = Xb; #SOLUÇÃO ÓTIMA
        disp('');
        disp('-------------------------------');
        disp('         SOLUÇÃO ÓTIMA:        '); 
        disp('-------------------------------');
        for j=1:n
           teste =[ "=> x", num2str(indiceVar(1,j)), ' = ',num2str(xSol(1, j))]; 
           disp(teste);  
        end
        disp('-------------------------------');
        custoOtimo = ['CUSTO OTIMO = ' num2str(z)];          
        disp(custoOtimo);    
        disp('-------------------------------'); 
        disp('');
        break;
    else    
        h = 0;  
        for i=1:m      
            if y(i,1)>0 #ELEMENTOS DA 'Direção Simplex' MAIORES DO QUE ZERO (POSITIVOS)
                h=h+1;
                auxTP = Xb(i,1)/y(i,1);#DIVISÃO DE ELEMENTOS DA 'Solução Básica'
                                       #POR ELEMENTOS DA 'Direção Simplex'
                                       
                TP(1, h) = auxTP;      #ARMAZENANDO RESULTADOS DA DIVISÃO ANTERIOS NO VETOR
                auxTP = [];          
            end
        end
        e = min(TP);      #TAMANHO DO PASSO
        
        x(indB) = Xb-y*e; #DETERMINA, ATRAVÉS DA 'Direção Simplex E 'Tamanho do Passo',
                          #VARIÁVEL BÁSICA QUE É ANULADA PARA A ESCOLHA DO 'índice da variável básica'
                          #QUE SAIRÁ DA BASE.
        #  x[indB] = [0 2 0 0 0]                
        
        indL= find(x(indB)==0, 1, "last"); #INDICE DO ELEMENTO A SAIR DA BASE
        saiDaBase = indB(1,indL);          #ELEMENTO QUE SAI DA PARTIÇÃO BÁSICA  
        
        indB(1, indL) = entraNaBase;#NOVA PARTIÇÃO BÁSICA
        indN(1, indK) =  saiDaBase; #NOVA PARTIÇÃO NÃO BÁSICA   
        
    end
    contAux = contAux+1;
endwhile
#FIM SIMPLEX
