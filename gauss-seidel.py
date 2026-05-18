from platform import system

import numpy as np

def P2R(abs: float, angle: float):
    angle = np.deg2rad(angle)  # Converter angulo de graus para radianos
    return abs*np.cos(angle) + 1j*abs*np.sin(angle)

def R2P(x: complex):
    return np.abs(x), np.angle(x)


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

        self.Zlines = {}  # Dicionário para armazenar as impedancias das linhas entre as barras
        for i in range(self.num_buses):
            for j in range(i+1, self.num_buses):
                self.Zlines[(i, j)] = 0.0 + 0.0j  # Impedancia da linha entre as barras i e j.

    def add_Zline(self, bus1: int, bus2: int, Z: complex) -> None:
        
        bus_index = (bus1, bus2)
        mirrowed_bus_index = (bus2, bus1)

        if bus_index in self.Zlines or mirrowed_bus_index in self.Zlines:
            self.Zlines[bus_index] = Z
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
            self.Qspec[bus] = Q
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

    def calculate_Ybus(self) -> None:

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


def gauss_seidel(sys: PowerSystem, tolerance: float=10e-4, max_iterations: int=100):
    V_iter = sys.V0.copy()  # Armazenar as tensoes anteriores para verificar a convergencia
    for iteration in range(max_iterations):
        for bus in range(sys.num_buses):
            if sys.bus_types[bus] == 'Slack':
                continue  # A barra slack tem tensao fixa, entao nao precisa ser atualizada
            if sys.bus_types[bus] == 'PQ':
                S_espec = sys.Pespec[bus] + 1j*sys.Qespec
                V_temp = S_espec / np.conj(V_iter[bus])  # Calcular a corrente injetada na barra
                for j in range(sys.num_buses):
                    if j != bus:
                        V_temp += sys.Ybus[bus, j] * V_iter[j]  # Somar as contribuicoes das outras barras
                V_iter[bus] = V_temp / sys.Ybus[bus, bus]  # Atualizar a tensao da barra usando a impedancia da barra
            elif sys.bus_types[bus] == 'PV':
                Q_calc = 
                


if __name__ == "__main__":
    pwsys = PowerSystem(['Slack', 'PQ', 'PV']) # onde o index das barras eh correspondente a sua posicao na lista, ou seja, barra 0 = PQ, barra 1 = PV e barra 2 = Slack
    pwsys.add_Zline(0, 1, 0.02 + 0.04j)
    pwsys.add_Zline(0, 2, 0.01 + 0.03j)
    pwsys.add_Zline(1, 2, 0.0125 + 0.025j)
    pwsys.calculate_Ybus()
    print(pwsys.Ybus)
    pwsys.add_bus_params(0, 1.0, 0.0)
    pwsys.add_bus_params(1, 1.0, 0.0, P=-4, Q=-2.5)
    pwsys.add_bus_params(2, 1.04, 0.0, P=2.0)
    gauss_seidel(pwsys)
    print("Tensões finais nas barras:")
    for bus in range(pwsys.num_buses):
        Vabs, Vangle = R2P(pwsys.V0[bus])
        print(f"Barra {bus}: |V| = {Vabs:.4f} p.u., Ângulo = {np.rad2deg(Vangle):.2f} graus")

