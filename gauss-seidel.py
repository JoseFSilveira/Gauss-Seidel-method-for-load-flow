import numpy as np

def P2R(abs: float, angle: float):
    angle = np.deg2rad(angle)  # Converter angulo de graus para radianos
    return abs*np.cos(angle) + 1j*abs*np.sin(angle)

def R2P(x: complex):
    return np.abs(x), np.rad2deg(np.angle(x))


class PowerSystem:

    def __init__(self, bus_types: list[str]) -> None:

        # Verificar se os tipos de barras informados sao validos
        possible_bus_types = ['PQ', 'PV', 'Slack']
        for bus_type in bus_types:
            if bus_type not in possible_bus_types:
                raise ValueError(f"Tipo de barra inválido: {bus_type}. Tipos válidos são: {possible_bus_types}")

        self.num_buses = len(bus_types)
        self.bus_types = bus_types
        self.Ybus = np.zeros((self.num_buses, self.num_buses), dtype=complex)
        self.V0 = np.ones(self.num_buses, dtype=complex)  # Magnitude das tensoes iniciais (1.0 p.u.)
        self.Pespec = np.zeros(self.num_buses, dtype=np.float64)  # Potencia especificada nas barras em pu
        self.Qespec = np.zeros(self.num_buses, dtype=np.float64)  # Potencia reativa especificada nas barras em pu
        self.delta = np.zeros(self.num_buses)  # Angulo das tensoes iniciais (0 degrees)

        self.Zlines = {}  # Dicionario para armazenar as impedancias das linhas entre as barras em p.u.
        self.B2_lines = {} # Dicionario para armazenar as susceptancias das linhas entre as barras em p.u.
        for i in range(self.num_buses):
            for j in range(i+1, self.num_buses):
                self.Zlines[(i, j)] = np.inf # Inicializar as impedancias das linhas com um valor infinito para indicar que as linhas ainda nao foram definidas
                self.B2_lines[(i, j)] = 0 # Inicializar as susceptancias das linhas com zero para indicar que as linhas ainda nao foram definidas

    def add_Zline(self, bus1: int, bus2: int, Z: complex, B2: float=0) -> None:
        
        if bus1 == bus2:
            raise ValueError("Não é possível adicionar uma linha entre a mesma barra. Verifique os índices das barras informados.")
        elif bus1 < 0 or bus1 >= self.num_buses or bus2 < 0 or bus2 >= self.num_buses:
            raise ValueError("Índices de barras inválidos. Verifique os índices das barras informados.")
        elif bus2 < bus1: # Garantir que bus1 sempre seja menor que bus2 para manter a consistencia do dicionario de impedancias das linhas
            bus1, bus2 = bus2, bus1

        bus_index = (bus1, bus2)
        mirrowed_bus_index = (bus2, bus1)

        if bus_index in self.Zlines or mirrowed_bus_index in self.Zlines:
            self.Zlines[bus_index] = Z
            self.B2_lines[bus_index] = B2
        else:
            raise ValueError(f"Não existe linha entre as barras {bus1} e {bus2}.")
        
    def add_bus_params(self, bus: int, Vabs: float, Vangle: float, P: float=None, Q: float=None) -> None:
        
        # Adicionar os parametros de tensao das barras
        self.V0[bus] = P2R(Vabs, Vangle)

        # Adicionar os parametros de potencia ativa e reativa das barras
        if self.bus_types[bus] == 'PQ':
            if P is None or Q is None:
                raise ValueError(f"Para barras do tipo PQ, é necessário especificar tanto a potência ativa (P) quanto a potência reativa (Q).")
            self.Pespec[bus] = P
            self.Qespec[bus] = Q
        elif self.bus_types[bus] == 'PV':
            if P is None:
                raise ValueError(f"Para barras do tipo PV, é necessário especificar a potência ativa (P).")
            self.Pespec[bus] = P
            self.Qespec[bus] = None  # Potencia reativa nao especificada para
        elif self.bus_types[bus] == 'Slack':
            self.Pespec[bus] = None  # Potencia ativa nao especificada para barra slack
            self.Qespec[bus] = None  # Potencia reativa nao especificada para barra slack
        else:
            raise ValueError(f"Tipo de barra inválido: {self.bus_types[bus]}. Tipos válidos são: ['PQ', 'PV', 'Slack']")

    def calculate_Ybus(self, verbose=False) -> None:

        # Zerar a matriz Ybus antes de calcular os valores, para evitar que valores antigos sejam mantidos caso a função seja chamada mais de uma vez
        self.Ybus = np.zeros((self.num_buses, self.num_buses), dtype=complex)
        
        # Verificar se as impedancias das linhas foram definidas corretamente
        for Zvalue in self.Zlines.values():
            if Zvalue == 0.0 + 0.0j:
                raise ValueError("Impedância de linha não definida corretamente. Verifique as linhas adicionadas.")

        for bus1 in range(self.num_buses):
            for bus2 in range(self.num_buses):
                if (bus1, bus2) in self.Zlines:
                    ybus = 1 / self.Zlines[(bus1, bus2)]
                    self.Ybus[bus1, bus2] = -ybus  # Elemento fora da diagonal e negativo
                    self.Ybus[bus2, bus1] = -ybus  # Elemento simetrico fora da diagonal
        for bus in range(self.num_buses):
            self.Ybus[bus, bus] = -np.sum(self.Ybus[bus, :])  # Elemento diagonal eh a soma dos elementos fora da diagonal, mas positivo
            self.Ybus[bus, bus] += np.sum([1j*B2 for (b1, b2), B2 in self.B2_lines.items() if bus in (b1, b2)])  # Adicionar as susceptancias shunt das linhas conectadas a barra
        
        if verbose:
            print(f"Matriz Ybus calculada:\n{self.Ybus}\n")


class GaussSeidel:
    def __init__(self, sys: PowerSystem) -> None:
        self.sys = sys
        self.V_iter = sys.V0.copy()  # Armazenar as tensoes anteriores para verificar a convergencia

    @staticmethod
    def calc_new_V(bus: int, V_iter: np.ndarray, Ybus: np.ndarray, S_espec: complex) -> complex:
        new_V = np.conj(S_espec) / np.conj(V_iter[bus])  # Calcular a corrente injetada na barra
        new_V -= np.sum(Ybus[bus, :] * V_iter, where=(np.arange(Ybus.shape[0]) != bus) )  # Somar as contribuicoes das outras barras
        return new_V / Ybus[bus, bus]  # Dividir pela admitancia da barra para obter a nova tensao


    def iterate(self, tolerance: float=1e-4, max_iterations: int=100, verbose: bool=False, coordinate_system: str='polar') -> None:
        '''
        Realiza as iteracoes do metodo de Gauss-Seidel para o sistema de potencia fornecido, atualizando as tensoes das barras a cada iteracao.
        Verbose habillita a impressao dos resultados de cada iteracao.
        coordinate_system define o sistema de coordenadas para a impressao dos resultados, podendo ser 'polar' ou 'rectangular'.
        '''

        if verbose and coordinate_system not in ['polar', 'rectangular']:
            raise ValueError(f"Sistema de coordenadas inválido: {coordinate_system}. Opções válidas são: ['polar', 'rectangular']")

        convergence_flag = False
        iteration = 0
        while not convergence_flag and iteration < max_iterations:
            if verbose:
                print((f"ITERAÇÃO {iteration+1}"))

            # Armazenar as tensoes anteriores para verificar a convergencia
            V_old = self.V_iter.copy()  # Armazenar as tensoes anteriores para verificar a convergencia

            for bus in range(self.sys.num_buses):
                
                if self.sys.bus_types[bus] == 'Slack':
                    pass  # A barra slack tem tensao fixa, entao nao precisa ser atualizada
                elif self.sys.bus_types[bus] == 'PQ':
                    S_espec = self.sys.Pespec[bus] + 1j*self.sys.Qespec[bus]  # Potencia especificada na barra
                    self.V_iter[bus] = self.calc_new_V(bus, self.V_iter, self.sys.Ybus, S_espec)
                elif self.sys.bus_types[bus] == 'PV':
                    Q_calc = -np.imag( np.conj(self.V_iter[bus]) * np.sum(self.sys.Ybus[bus, :] * self.V_iter) )  # Calcular a potencia reativa injetada na barra
                    S_espec = self.sys.Pespec[bus] + 1j*Q_calc  # Potencia especificada na barra, usando a potencia ativa especificada e a potencia reativa calculada
                    V_calc = self.calc_new_V(bus, self.V_iter, self.sys.Ybus, S_espec)
                    _, Vangle_calc = R2P(V_calc)
                    self.V_iter[bus] = P2R(np.abs(self.sys.V0[bus]), Vangle_calc)  # Manter a magnitude da tensao especificada e atualizar o angulo calculado
                else:
                    raise ValueError(f"Tipo de barra inválido: {self.sys.bus_types[bus]}. Tipos válidos são: ['PQ', 'PV', 'Slack']")
                
                # Escrever o valor da tensao calculada para a barra
                if verbose:
                    if coordinate_system == 'rectangular':
                        print(f"Barra {bus} ({self.sys.bus_types[bus]}): Tensão calculada = {self.V_iter[bus]} p.u.")
                    elif coordinate_system == 'polar':
                        Vabs, Vangle = R2P(self.V_iter[bus])
                        print(f"Barra {bus} ({self.sys.bus_types[bus]}): Tensão calculada = {Vabs} | {Vangle}° p.u.")

            # Verificar a convergencia comparando as tensoes atuais com as tensoes anteriores
            max_diff = np.max(np.abs(self.V_iter - V_old))
            if max_diff < tolerance:
                print(f"Convergencia atingida após {iteration+1} iterações.")
                convergence_flag = True
            iteration += 1 # aumentar o contador de iteracoes
            if iteration >= max_iterations and not convergence_flag:
                print(f"Limite máximo de iterações atingido sem convergência. Diferença máxima: {max_diff}")
            print()  # Linha em branco para separar as iteracoes


if __name__ == "__main__":

    # Criar o sistema de potencia
    pwsys = PowerSystem(['Slack', 'PV', 'PQ', 'PQ']) # onde o index das barras eh correspondente a sua posicao na lista, ou seja, barra 0 = PQ, barra 1 = PV e barra 2 = Slack

    # Adicionando as impedancias das linhas entre as barras
    pwsys.add_Zline(0, 1, 0.15 + 0.4j, B2=0.04)
    pwsys.add_Zline(0, 2, 0.1 + 0.3j, B2=0.05)
    pwsys.add_Zline(0, 3, 0.15 + 0.6j, B2=0.04)
    pwsys.add_Zline(1, 2, 0.07 + 0.25j, B2=0.03)
    pwsys.add_Zline(2, 3, 0.09 + 0.3j, B2=0.04)

    # Calcular a matriz Ybus
    print() # Adicionar espaco para melhor visualizacao dos resultados
    pwsys.calculate_Ybus(verbose=True)

    # Adicionar os parametros das barras (tensao, potencia ativa e potencia reativa)
    pwsys.add_bus_params(0, 1.06, 0.0)
    pwsys.add_bus_params(1, 1.02, 0.0, P=1.2)
    pwsys.add_bus_params(2, 1.0, 0.0, P=-0.7, Q=-0.5)
    pwsys.add_bus_params(3, 1.0, 0.0, P=-0.6, Q=-0.3)

    # Calcular as tensoes das barras usando o metodo de Gauss-Seidel
    gs = GaussSeidel(pwsys)
    gs.iterate(tolerance=1e-4, max_iterations=100, verbose=True, coordinate_system='rectangular')