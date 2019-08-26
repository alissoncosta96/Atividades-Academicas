#DUPLA: FRANCISCO ALISSON E FRANCISCO REGINALDO


############## ARMAZENAMENTO DE DADOS - INICIO ################
C = input('Digite a Fun��o Objetivo na forma vetorial [X1 X2 ... XN]: ');#ENTRADA 'fun��o objetivo'
disp('');
m = input('Digite a quantidade de restri��es: '); #ENTRADA 'quantidade restri��es'

for i= 1:m
    Ab(i,:)=[input('Digite a restri��o na forma vetorial: ')]; #ARMAZENA AS restri��es JUNTAMENTE 
                                                               #COM OS 'termos ap�s a igualdade'
end

nAb = size(Ab, 2); #QUANTIDADE DE COLUNAS DE Ab

for i=1:nAb-1   
     A(:,i) = Ab(:,i); #RETIRA A ultima coluna 'b' 
                       #DA MATRIZ Ab SEPARANDO 'A' e 'b' - MATRIZ 'A'
end

b = Ab(:, nAb); #ARMAZENA A MATRIZ 'b' 

n = size(A,2); #QUANTIDADE DE COLUNAS DE 'A'
varIndep = n - m; #QUANTIDADE DE VARI�VEIS INDEPENDENTES QUE PODEM SER FIXADAS COMO ZERO

k=0;    
for j=1:n 
    k = k+1;
    indiceVar(1, j) = k;#ATRIBUI '�ndices' A TODAS AS 'colunas (vari�veis)'
end  

Xb = zeros(1,n-varIndep);#INICIALIZANDO AS 'vari�veis b�sicas' COM 'zero'
############# ARMAZENAMENTO DE DADOS - FIM #################
 

#IN�CIO DO SIMPLEX
g = 1;
contAux = 1; 
while g==1
  p = 1;
    while p==1   
        if contAux==1#VERIFICA SE os �ndices B e N 
            if isempty(find(Xb<0))# VERICA SE TODOS OS ELEMENTOS 
                                  # DE 'Xb' S�O POSITIVOS (FACTIBILIDADE)
                                  
                indN = indiceVar(1:varIndep);#DEFINE OS '�ndices das vari�veis b�sicas'
                                             #PARTI��O N�O B�SICA
                                             
                indB = indiceVar(varIndep+1:n);#DEFINE OS '�ndices das vari�veis n�o b�sicas'
                                               #PARTI��O B�SICA
                
            else #CASO Xb 'n�o seja fact�vel'
                indN = indiceVar(n-varIndep+1:n); #DEFINE 'novos �ndices das vari�veis b�sicas'
                indB = indiceVar(1:n-varIndep);   #DEFINE 'novos �ndices das vari�veis n�o b�sicas'                
            end  
        end 
        
        B = A(:, indB); #MATRIZ 'B'
        N = A(:, indN); #MATRIZ 'N'
            
        
       Xb=inv(B)*b;   
       if isempty(find(Xb<0)) #VERIFICA A 'factibilidade' DA SOLU��O B�SICA     
           p=0;
       end    
    endwhile

    Cb = C(indB);#'fun��o objetivo' ESCRITA NOS    
                 #'�ndices das vari�veis b�sicas'
                 
    Cn = C(indN);#'fun��o objetivo' ESCRITA NOS 
                 #'�ndices das vari�veis n�o b�sicas'

    mSimplex =  Cb*inv(B); #'Multiplicador Simplex': NECESS�RIO PARA
                           # CALCULAR OS 'Custos Reduzidos'
    
    
    q=1;   
    
    while q<=numel(indN)
        auxRed = Cn(1,q) - ((mSimplex)*(N(:, q))); #  CALCULA OS 'Custos Reduzidos'
        cRed(1,q) = auxRed; # ARMAZENA OS 'Custos Reduzidos' NO VETOR
        auxRed = [];
        q++;     
    endwhile
    
    [menorCR, indK] = min(cRed);# ARMAZENA O 'menor custo reduzido' COM SEU '�ndice'
                                # DETERMINA��O DO �NDICE A ENTRAR NA BASE  
                                
    entraNaBase =  indN(1,indK);# �NDICE QUE VAI ENTRAR NA 'base'
    
    y = inv(B)*N(:,indK);#  DIRE��O SIMPLEX
    
    if isempty(find(y>0))#  VERIFICA SE TODOS O ELEMENTOS 
                         #  DA 'Dire��o Simplex' S�O NEGATIVOS (SOLU��O ILIMITADA)
          disp('-------------------------------------------------------------------------');
          disp("A SOLU��O �TIMA � ILIMITADA COM A FUN��O CUSTO AUMENTANDO PARA O INFINITO");
          disp('-------------------------------------------------------------------------');
          break;
    end
    
    if isempty(find(cRed<0))# VERIFICA SE TODOS OS ELEMENTOS
                            # DO 'Custo Reduzido' S�O POSITIVOS (SOLU��O �TIMA) 
        g=0;
        z = Cb*Xb; # CUSTO �TIMO
        xSol = zeros(1, n);
        xSol(indB) = Xb; #SOLU��O �TIMA
        disp('');
        disp('-------------------------------');
        disp('         SOLU��O �TIMA:        '); 
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
            if y(i,1)>0 #ELEMENTOS DA 'Dire��o Simplex' MAIORES DO QUE ZERO (POSITIVOS)
                h=h+1;
                auxTP = Xb(i,1)/y(i,1);#DIVIS�O DE ELEMENTOS DA 'Solu��o B�sica'
                                       #POR ELEMENTOS DA 'Dire��o Simplex'
                                       
                TP(1, h) = auxTP;      #ARMAZENANDO RESULTADOS DA DIVIS�O ANTERIOS NO VETOR
                auxTP = [];          
            end
        end
        e = min(TP);      #TAMANHO DO PASSO
        
        x(indB) = Xb-y*e; #DETERMINA, ATRAV�S DA 'Dire��o Simplex E 'Tamanho do Passo',
                          #VARI�VEL B�SICA QUE � ANULADA PARA A ESCOLHA DO '�ndice da vari�vel b�sica'
                          #QUE SAIR� DA BASE.
        #  x[indB] = [0 2 0 0 0]                
        
        indL= find(x(indB)==0, 1, "last"); #INDICE DO ELEMENTO A SAIR DA BASE
        saiDaBase = indB(1,indL);          #ELEMENTO QUE SAI DA PARTI��O B�SICA  
        
        indB(1, indL) = entraNaBase;#NOVA PARTI��O B�SICA
        indN(1, indK) =  saiDaBase; #NOVA PARTI��O N�O B�SICA   
        
    end
    contAux = contAux+1;
endwhile
#FIM SIMPLEX