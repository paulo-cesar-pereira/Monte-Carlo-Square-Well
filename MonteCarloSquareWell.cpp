#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <cmath>
using namespace std;

	const int N = 64;								// Número de partículas
	const int nx = sqrt(N);							// Número de linhas/colunas da rede (para rede quadrada)
	const float a0 = 1.5f;							// Parâmetro de rede
	const float L = nx*a0;							// Tamanho da caixa
	const float mL = L*0.5f;						// Metade do tamanho da caixa
	const float alpha = 1.65f;						// Alcance do potencial de interação
	float drmax = a0*0.05f;							// Passo máximo de deslocamento
	const float auxangulo = 2.0f * M_PI / RAND_MAX;	// Float auxiliar para cálculo do ângulo aleatório	
	float auxdr = drmax / RAND_MAX;					// Float auxiliar para cálculo do deslocamento aleatório
	float T = 10.0f;								// Temperatura
	int maxpassos = 1000000;						// Número máximo de passos Monte Carlo


// Estrutura para representar um vetor 2D
struct Float2 {
    float x;
    float y;

    // Construtor
    Float2(float x = 0.0f, float y = 0.0f) : x(x), y(y) {}

    // Operador de soma
    Float2 operator+(const Float2& other) const {
        return {x + other.x, y + other.y};
    }
    
    // Operador de subtração
    Float2 operator-(const Float2& other) const {
        return {x - other.x, y - other.y};
    }
    
    // Módulo
    float modulo(){
        return sqrtf(x * x + y * y);
    }
};

// Função que corrige a posição, caso ela caia fora da caixa
void corrigePosicao(Float2& pos){
    if(pos.x > L) pos.x -= L;
    else if(pos.x < 0.0f) pos.x += L;
    if(pos.y > L) pos.y -= L;
    else if(pos.y < 0.0f) pos.y += L;
}

// Função que corrige a distância entre duas partículas, no caso de haver partícula imagem mais próxima
void corrigeDistancia(Float2& dist){
    if(dist.x > mL) dist.x -= L;
    else if(dist.x < -mL) dist.x += L;
    if(dist.y > mL) dist.y -= L;
    else if(dist.y < -mL) dist.y += L;
}

// Função que cria um delocamento aleatório
Float2 deslocamentoAleatorio(){
	float angulo = rand() * auxangulo;
	float dr = rand() * auxdr;
	return { dr * cos(angulo), dr * sin(angulo) };
}


// Estrutura que gerencia um conjunto de pontos com alocação dinâmica
struct Posicoes {
    int tamanho;
    Float2* pontos;
    int* energias;
    Float2 novapos;
	int novaenergia;

    // Construtor: aloca memória
    Posicoes(int n) : tamanho(n) {
        pontos = new Float2[n];
        energias = new int[n];
        int nx = sqrt(n);
    	for(int i=0;i<nx;i++){
    		for(int j=0;j<nx;j++){
    			pontos[i*nx+j].x = (i+0.5f)*a0;
    			pontos[i*nx+j].y = (j+0.5f)*a0;
    			energias[i*nx+j] = 0;
    		}
    	}
    	novapos = {0.0f, 0.0f};
		novaenergia = 0;
    	for (int i = 0; i < tamanho; i++) energias[i] = variacaoEnergia(i, {0.0f, 0.0f});
    }

    // Destrutor: libera memória
    ~Posicoes() {
        delete[] pontos;
        delete[] energias;
    }
    
    // Calcula a variação de energia de interação de uma dada partícula com as outras, considerando um deslocamento em sua posição
    int variacaoEnergia(int i, const Float2& desl) {
        novapos = pontos[i] + desl;
        corrigePosicao(novapos);
        Float2 distanciaij;
        float distancia;
        novaenergia = 0;
        for (int j = 0; j < tamanho; j++) {
            if(j != i){
                distanciaij = novapos - pontos[j];
                corrigeDistancia(distanciaij);
                distancia = distanciaij.modulo();
                if(distancia >= 1.0f && distancia < alpha) novaenergia -= 1;
            }
        }
        return novaenergia - energias[i];
    }
    // Atualiza a posição de uma partícula
	void atualizaPosicao(int i){
		pontos[i] = novapos;
		energias[i] = novaenergia;
	}
};


int main (){

	// Inicializa o sistema
    Posicoes* particulas = new Posicoes(N);

	// Exporta as posições iniciais
	ofstream posiniciais("posiniciais.dat");
	for(int i=0;i<N;i++){
		posiniciais << particulas->pontos[i].x << "   " << particulas->pontos[i].y << "   " << particulas->energias[i] << endl;
	}
	posiniciais.close();

	// Escolhe o seed
	srand(time(NULL));
	
	// Declara variáveis úteis
	int ia, varenergia;
	float prob, sorteio;

	// Passos de Monte Carlo
	for(int p = 0; p < maxpassos; p++){
		
		// Sorteia uma partícula
		ia = rand() % N;

		// Calcula a variação de energia para um deslocamento aleatório dessa partícula
		varenergia = particulas->variacaoEnergia(ia, deslocamentoAleatorio(drmax));
		
		// Algoritmo de Metropolis aceita a nova configuração, caso a energia não tenha aumentado
		if(varenergia <= 0) particulas->atualizaPosicao(ia);
		else{
			// Caso a energia tenha aumentado, a probabilidade da nova configuração ser aceita depende da temperatura
			sorteio = static_cast<float>(rand())/RAND_MAX;
			prob = expf(-varenergia / T);
			if(sorteio <= prob) particulas->atualizaPosicao(ia);
		}
	}


	// Exporta as posições finais:
	ofstream posfinais("posfinais.dat");
	for(int i=0;i<N;i++){
		posfinais << particulas->pontos[i].x << "   " << particulas->pontos[i].y << "   " << particulas->energias[i] << endl;
	}
	posfinais.close();

	delete particulas;
	return 0;

}
