import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d import Axes3D
import json

# Tentar importar temas modernos para o Tkinter, caso contrário, usar o padrão
try:
    from ttkthemes import ThemedTk
    TEMA_DISPONIVEL = True
except ImportError:
    TEMA_DISPONIVEL = False

class AplicativoVisualizadorGeometriaAvancada:
    def __init__(self, raiz):
        self.raiz = raiz
        self.raiz.title("Visualizador Avançado de Geometria Analítica")
        self.raiz.state('zoomed')  # Maximizar janela

        # Criar notebook para diferentes visualizações
        self.caderno = ttk.Notebook(raiz)
        self.caderno.pack(fill=tk.BOTH, expand=True)

        # Páginas principais
        self.pagina_geometria = ttk.Frame(self.caderno)
        self.pagina_mapa_conceitual = ttk.Frame(self.caderno)

        self.caderno.add(self.pagina_geometria, text="Visualização Geométrica")
        self.caderno.add(self.pagina_mapa_conceitual, text="Mapa Conceitual")

        # Configurar páginas
        self.configurar_pagina_geometria()
        self.configurar_pagina_mapa_conceitual()

        # Estado de interação para visualização 3D
        self.arrastando = False
        self.ultimo_x = 0
        self.ultimo_y = 0
        self.rotacao = [30, 45, 0]
        self.translacao = [0, 0, 0]

    def configurar_pagina_geometria(self):
        # Frame principal para geometria
        self.frame_geometria = ttk.PanedWindow(self.pagina_geometria, orient=tk.HORIZONTAL)
        self.frame_geometria.pack(fill=tk.BOTH, expand=True)

        # Painel de controle
        self.frame_controle = ttk.Frame(self.frame_geometria)
        self.frame_geometria.add(self.frame_controle, weight=1)

        # Área de visualização
        self.frame_visualizacao = ttk.Frame(self.frame_geometria)
        self.frame_geometria.add(self.frame_visualizacao, weight=3)

        # Configurar controles e plotagem
        self.configurar_controles()
        self.configurar_controles_parametros()
        self.configurar_plotagem()

    def configurar_controles(self):
        # Seleção de objeto geométrico
        ttk.Label(self.frame_controle, text="Objeto Geométrico:").pack(pady=5)
        self.objeto_var = tk.StringVar(value="Parábola")
        objetos = [
            "Parábola", "Elipse", "Hipérbole", "Elipsóide",
            "Hiperbolóide (1 folha)", "Hiperbolóide (2 folhas)",
            "Parabolóide Hiperbólico", "Parabolóide Elíptico",
            "Superfície Cônica", "Superfície Cilíndrica", "Superfície Esférica"
        ]

        self.menu_objeto = ttk.Combobox(
            self.frame_controle,
            textvariable=self.objeto_var,
            values=objetos,
            state="readonly"
        )
        self.menu_objeto.pack(pady=5)
        self.menu_objeto.bind('<<ComboboxSelected>>', self.atualizar_plotagem)

    def configurar_controles_parametros(self):
        # Frame para parâmetros
        self.frame_parametros = ttk.LabelFrame(self.frame_controle, text="Parâmetros")
        self.frame_parametros.pack(fill=tk.X, padx=5, pady=5)

        # Parâmetros específicos para cada objeto
        self.parametros_widgets = {}

        # Definição de parâmetros para todos os objetos geométricos
        self.parametros_definidos = {
            'parábola': {
                'a': {'valor': 1.0, 'intervalo': (-5, 5), 'rotulo': 'Abertura'},
                'h': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Deslocamento X'},
                'k': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Deslocamento Y'}
            },
            'elipse': {
                'a': {'valor': 2.0, 'intervalo': (0.1, 5), 'rotulo': 'Semieixo Maior'},
                'b': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Semieixo Menor'},
                'h': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Centro X'},
                'k': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Centro Y'},
                'rotacao': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação (°)'}
            },
            'hipérbole': {
                'a': {'valor': 2.0, 'intervalo': (0.1, 5), 'rotulo': 'Semieixo Transversal'},
                'b': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Semieixo Conjugado'},
                'h': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Centro X'},
                'k': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Centro Y'},
                'rotacao': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação (°)'}
            },
            'elipsóide': {
                'a': {'valor': 2.0, 'intervalo': (0.1, 5), 'rotulo': 'Semieixo X'},
                'b': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Semieixo Y'},
                'c': {'valor': 1.5, 'intervalo': (0.1, 5), 'rotulo': 'Semieixo Z'},
                'h': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Centro X'},
                'k': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Centro Y'},
                'l': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Centro Z'},
                'rotacao_x': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação X (°)'},
                'rotacao_y': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação Y (°)'},
                'rotacao_z': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação Z (°)'}
            },
            'hiperbolóide (1 folha)': {
                'a': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Semieixo X'},
                'b': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Semieixo Y'},
                'c': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Semieixo Z'},
                'h': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Centro X'},
                'k': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Centro Y'},
                'l': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Centro Z'},
                'rotacao_x': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação X (°)'},
                'rotacao_y': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação Y (°)'},
                'rotacao_z': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação Z (°)'}
            },
            'hiperbolóide (2 folhas)': {
                'a': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Semieixo X'},
                'b': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Semieixo Y'},
                'c': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Semieixo Z'},
                'h': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Centro X'},
                'k': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Centro Y'},
                'l': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Centro Z'},
                'rotacao_x': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação X (°)'},
                'rotacao_y': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação Y (°)'},
                'rotacao_z': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação Z (°)'}
            },
            'parabolóide hiperbólico': {
                'a': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Coeficiente X²'},
                'b': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Coeficiente Y²'},
                'h': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Deslocamento X'},
                'k': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Deslocamento Y'},
                'l': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Deslocamento Z'}
            },
            'parabolóide elíptico': {
                'a': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Coeficiente X²'},
                'b': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Coeficiente Y²'},
                'h': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Deslocamento X'},
                'k': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Deslocamento Y'},
                'l': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Deslocamento Z'}
            },
            'superfície cônica': {
                'a': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Coeficiente X'},
                'b': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Coeficiente Y'},
                'c': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Coeficiente Z'},
                'h': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Deslocamento X'},
                'k': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Deslocamento Y'},
                'l': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Deslocamento Z'},
                'rotacao_x': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação X (°)'},
                'rotacao_y': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação Y (°)'},
                'rotacao_z': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação Z (°)'}
            },
            'superfície cilíndrica': {
                'a': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Raio'},
                'h': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Deslocamento X'},
                'k': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Deslocamento Y'},
                'l': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Deslocamento Z'},
                'rotacao_x': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação X (°)'},
                'rotacao_y': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação Y (°)'},
                'rotacao_z': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação Z (°)'}
            },
            'superfície esférica': {
                'a': {'valor': 1.0, 'intervalo': (0.1, 5), 'rotulo': 'Raio'},
                'h': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Centro X'},
                'k': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Centro Y'},
                'l': {'valor': 0.0, 'intervalo': (-5, 5), 'rotulo': 'Centro Z'},
                'rotacao_x': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação X (°)'},
                'rotacao_y': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação Y (°)'},
                'rotacao_z': {'valor': 0.0, 'intervalo': (-180, 180), 'rotulo': 'Rotação Z (°)'}
            }
        }

        self.criar_widgets_parametros()

    def criar_widgets_parametros(self):
        objeto_atual = self.objeto_var.get().lower()

        # Limpar widgets existentes
        for widgets in self.parametros_widgets.values():
            for widget in widgets:
                widget.destroy()
        self.parametros_widgets.clear()

        # Selecionar parâmetros com base no objeto
        parametros = self.parametros_definidos.get(objeto_atual, {})

        # Criar novos widgets baseados no objeto selecionado
        for nome_param, info_param in parametros.items():
            frame = ttk.Frame(self.frame_parametros)
            frame.pack(fill=tk.X, padx=5, pady=2)

            rotulo = ttk.Label(frame, text=info_param['rotulo'])
            rotulo.pack(side=tk.LEFT)

            var = tk.DoubleVar(value=info_param['valor'])
            entrada = ttk.Entry(frame, textvariable=var, width=8)
            entrada.pack(side=tk.LEFT, padx=5)

            escala = ttk.Scale(
                frame,
                from_=info_param['intervalo'][0],
                to=info_param['intervalo'][1],
                variable=var,
                orient=tk.HORIZONTAL,
                command=lambda event, nome=nome_param: self.atualizar_plotagem()
            )
            escala.pack(side=tk.LEFT, fill=tk.X, expand=True)

            self.parametros_widgets[nome_param] = [frame, rotulo, entrada, escala]
            var.trace_add('write', lambda *args, nome=nome_param: self.atualizar_plotagem())

    def configurar_plotagem(self):
        # Criar figura do Matplotlib
        self.figura = Figure(figsize=(8, 6), dpi=100)
        self.canvas_mpl = FigureCanvasTkAgg(self.figura, master=self.frame_visualizacao)
        self.canvas_mpl.draw()
        self.canvas_mpl.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Barra de ferramentas de navegação
        from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
        self.toolbar = NavigationToolbar2Tk(self.canvas_mpl, self.frame_visualizacao)
        self.toolbar.update()

        # Conectar eventos do mouse para interação 3D
        self.canvas_mpl.mpl_connect('button_press_event', self.on_mouse_press)
        self.canvas_mpl.mpl_connect('button_release_event', self.on_mouse_release)
        self.canvas_mpl.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas_mpl.mpl_connect('scroll_event', self.on_scroll)

        # Plotagem inicial
        self.atualizar_plotagem()

    def atualizar_plotagem(self, evento=None):
        self.figura.clear()

        # Obter parâmetros do objeto atual
        objeto_atual = self.objeto_var.get().lower()
        parametros = self.obter_parametros_atual()

        # Determinar se a plotagem será 2D ou 3D
        if any(x in objeto_atual for x in [
            "parabolóide", "hiperbolóide", "superfície", "elipsóide"
        ]):
            ax = self.figura.add_subplot(111, projection='3d')
            ax.view_init(elev=self.rotacao[0], azim=self.rotacao[1])
        else:
            ax = self.figura.add_subplot(111)

        try:
            # Plotar objeto selecionado com os parâmetros atuais
            if objeto_atual == "parábola":
                self.plotar_parabola(ax, parametros)
            elif objeto_atual == "elipse":
                self.plotar_elipse(ax, parametros)
            elif objeto_atual == "hipérbole":
                self.plotar_hiperbole(ax, parametros)
            elif objeto_atual == "elipsóide":
                self.plotar_elipsóide(ax, parametros)
            elif "hiperbolóide (1 folha)" in objeto_atual:
                self.plotar_hiperboloide_uma_folha(ax, parametros)
            elif "hiperbolóide (2 folhas)" in objeto_atual:
                self.plotar_hiperboloide_duas_folhas(ax, parametros)
            elif "parabolóide hiperbólico" in objeto_atual:
                self.plotar_paraboloide_hiperbolico(ax, parametros)
            elif "parabolóide elíptico" in objeto_atual:
                self.plotar_paraboloide_eliptico(ax, parametros)
            elif "superfície cônica" in objeto_atual:
                self.plotar_superficie_conica(ax, parametros)
            elif "superfície cilíndrica" in objeto_atual:
                self.plotar_superficie_cilindrica(ax, parametros)
            elif "superfície esférica" in objeto_atual:
                self.plotar_superficie_esferica(ax, parametros)

            # Configurações comuns do plot
            if not isinstance(ax, Axes3D):
                ax.grid(True)
                ax.set_aspect('equal')
                ax.set_xlim(-10 + self.translacao[0], 10 + self.translacao[0])
                ax.set_ylim(-10 + self.translacao[1], 10 + self.translacao[1])
            else:
                ax.grid(True)
                ax.set_box_aspect([1, 1, 1])

            self.figura.tight_layout()
            self.canvas_mpl.draw()

        except Exception as e:
            messagebox.showerror("Erro", f"Erro ao plotar: {str(e)}")

    def obter_parametros_atual(self):
        """Obtém os parâmetros atuais do objeto selecionado"""
        parametros = {}
        for nome_param, widgets in self.parametros_widgets.items():
            try:
                valor = float(widgets[2].get())  # Pega o valor do Entry
                parametros[nome_param] = valor
            except (ValueError, IndexError):
                parametros[nome_param] = 1.0  # Valor padrão em caso de erro
        return parametros

    # Funções para plotar diferentes objetos geométricos
    def plotar_parabola(self, ax, params):
        a = params.get('a', 1.0)
        h = params.get('h', 0.0)
        k = params.get('k', 0.0)

        x = np.linspace(-10, 10, 400)
        y = a * (x - h) ** 2 + k

        ax.plot(x, y, 'b-', linewidth=2, label='Parábola')
        ax.plot([h], [k], 'ro', label='Vértice')  # Marca o vértice
        ax.legend()

    def plotar_elipse(self, ax, params):
        a = params.get('a', 2.0)
        b = params.get('b', 1.0)
        h = params.get('h', 0.0)
        k = params.get('k', 0.0)
        rotacao = params.get('rotacao', 0.0)

        t = np.linspace(0, 2 * np.pi, 400)
        x = a * np.cos(t)
        y = b * np.sin(t)

        # Aplicar rotação
        theta = np.radians(rotacao)
        x_rot = x * np.cos(theta) - y * np.sin(theta)
        y_rot = x * np.sin(theta) + y * np.cos(theta)

        # Aplicar translação
        x_final = x_rot + h
        y_final = y_rot + k

        ax.plot(x_final, y_final, 'g-', linewidth=2, label='Elipse')
        ax.plot([h], [k], 'ro', label='Centro')
        ax.legend()

    def plotar_hiperbole(self, ax, params):
        a = params.get('a', 2.0)
        b = params.get('b', 1.0)
        h = params.get('h', 0.0)
        k = params.get('k', 0.0)
        rotacao = params.get('rotacao', 0.0)

        t = np.linspace(-10, 10, 400)
        x = a * np.cosh(t)
        y = b * np.sinh(t)

        # Aplicar rotação
        theta = np.radians(rotacao)
        x_rot = x * np.cos(theta) - y * np.sin(theta)
        y_rot = x * np.sin(theta) + y * np.cos(theta)

        # Aplicar translação
        x_final = x_rot + h
        y_final = y_rot + k

        ax.plot(x_final, y_final, 'm-', linewidth=2, label='Ramos da Hiperbole')
        ax.plot(-x_final, y_final, 'm-', linewidth=2)
        ax.plot([h], [k], 'ro', label='Centro')
        ax.legend()

    def plotar_elipsóide(self, ax, params):
        a = params.get('a', 2.0)
        b = params.get('b', 1.0)
        c = params.get('c', 1.5)
        h = params.get('h', 0.0)
        k = params.get('k', 0.0)
        l = params.get('l', 0.0)
        rotacao_x = params.get('rotacao_x', 0.0)
        rotacao_y = params.get('rotacao_y', 0.0)
        rotacao_z = params.get('rotacao_z', 0.0)

        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        U, V = np.meshgrid(u, v)
        X = a * np.outer(np.cos(U), np.sin(V)) + h
        Y = b * np.outer(np.sin(U), np.sin(V)) + k
        Z = c * np.outer(np.ones_like(U), np.cos(V)) + l

        # Aplicar rotações (opcional)
        if rotacao_x != 0 or rotacao_y != 0 or rotacao_z != 0:
            # Rotações em ordem Z, Y, X
            rot_z = np.radians(rotacao_z)
            rot_y = np.radians(rotacao_y)
            rot_x = np.radians(rotacao_x)

            # Matriz de rotação Z
            Rz = np.array([
                [np.cos(rot_z), -np.sin(rot_z), 0],
                [np.sin(rot_z),  np.cos(rot_z), 0],
                [0,             0,             1]
            ])

            # Matriz de rotação Y
            Ry = np.array([
                [ np.cos(rot_y), 0, np.sin(rot_y)],
                [0,              1, 0],
                [-np.sin(rot_y), 0, np.cos(rot_y)]
            ])

            # Matriz de rotação X
            Rx = np.array([
                [1, 0,              0],
                [0, np.cos(rot_x), -np.sin(rot_x)],
                [0, np.sin(rot_x),  np.cos(rot_x)]
            ])

            # Aplicar rotações
            for i in range(X.shape[0]):
                for j in range(X.shape[1]):
                    vec = np.array([X[i, j], Y[i, j], Z[i, j]])
                    vec_rot = Rz @ Ry @ Rx @ vec
                    X[i, j], Y[i, j], Z[i, j] = vec_rot

        ax.plot_surface(X, Y, Z, color='cyan', alpha=0.6, edgecolor='none', label='Elipsóide')

    def plotar_hiperboloide_uma_folha(self, ax, params):
        a = params.get('a', 1.0)
        b = params.get('b', 1.0)
        c = params.get('c', 1.0)
        h = params.get('h', 0.0)
        k = params.get('k', 0.0)
        l = params.get('l', 0.0)
        rotacao_x = params.get('rotacao_x', 0.0)
        rotacao_y = params.get('rotacao_y', 0.0)
        rotacao_z = params.get('rotacao_z', 0.0)

        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(-2, 2, 100)
        U, V = np.meshgrid(u, v)
        X = a * np.cosh(V) * np.cos(U) + h
        Y = b * np.cosh(V) * np.sin(U) + k
        Z = c * np.sinh(V) + l

        # Aplicar rotações (opcional)
        if rotacao_x != 0 or rotacao_y != 0 or rotacao_z != 0:
            # Rotações em ordem Z, Y, X
            rot_z = np.radians(rotacao_z)
            rot_y = np.radians(rotacao_y)
            rot_x = np.radians(rotacao_x)

            # Matriz de rotação Z
            Rz = np.array([
                [np.cos(rot_z), -np.sin(rot_z), 0],
                [np.sin(rot_z),  np.cos(rot_z), 0],
                [0,             0,             1]
            ])

            # Matriz de rotação Y
            Ry = np.array([
                [ np.cos(rot_y), 0, np.sin(rot_y)],
                [0,              1, 0],
                [-np.sin(rot_y), 0, np.cos(rot_y)]
            ])

            # Matriz de rotação X
            Rx = np.array([
                [1, 0,              0],
                [0, np.cos(rot_x), -np.sin(rot_x)],
                [0, np.sin(rot_x),  np.cos(rot_x)]
            ])

            # Aplicar rotações
            for i in range(X.shape[0]):
                for j in range(X.shape[1]):
                    vec = np.array([X[i, j], Y[i, j], Z[i, j]])
                    vec_rot = Rz @ Ry @ Rx @ vec
                    X[i, j], Y[i, j], Z[i, j] = vec_rot

        ax.plot_surface(X, Y, Z, color='yellow', alpha=0.6, edgecolor='none', label='Hiperbolóide (1 folha)')

    def plotar_hiperboloide_duas_folhas(self, ax, params):
        a = params.get('a', 1.0)
        b = params.get('b', 1.0)
        c = params.get('c', 1.0)
        h = params.get('h', 0.0)
        k = params.get('k', 0.0)
        l = params.get('l', 0.0)
        rotacao_x = params.get('rotacao_x', 0.0)
        rotacao_y = params.get('rotacao_y', 0.0)
        rotacao_z = params.get('rotacao_z', 0.0)

        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(-2, 2, 100)
        U, V = np.meshgrid(u, v)
        X = a * np.sinh(V) * np.cos(U) + h
        Y = b * np.sinh(V) * np.sin(U) + k
        Z = c * np.cosh(V) + l

        # Aplicar rotações (opcional)
        if rotacao_x != 0 or rotacao_y != 0 or rotacao_z != 0:
            # Rotações em ordem Z, Y, X
            rot_z = np.radians(rotacao_z)
            rot_y = np.radians(rotacao_y)
            rot_x = np.radians(rotacao_x)

            # Matriz de rotação Z
            Rz = np.array([
                [np.cos(rot_z), -np.sin(rot_z), 0],
                [np.sin(rot_z),  np.cos(rot_z), 0],
                [0,             0,             1]
            ])

            # Matriz de rotação Y
            Ry = np.array([
                [ np.cos(rot_y), 0, np.sin(rot_y)],
                [0,              1, 0],
                [-np.sin(rot_y), 0, np.cos(rot_y)]
            ])

            # Matriz de rotação X
            Rx = np.array([
                [1, 0,              0],
                [0, np.cos(rot_x), -np.sin(rot_x)],
                [0, np.sin(rot_x),  np.cos(rot_x)]
            ])

            # Aplicar rotações
            for i in range(X.shape[0]):
                for j in range(X.shape[1]):
                    vec = np.array([X[i, j], Y[i, j], Z[i, j]])
                    vec_rot = Rz @ Ry @ Rx @ vec
                    X[i, j], Y[i, j], Z[i, j] = vec_rot

        ax.plot_surface(X, Y, Z, color='orange', alpha=0.6, edgecolor='none', label='Hiperbolóide (2 folhas)')

        # Plotar a segunda folha
        X2 = -X
        Y2 = Y
        Z2 = -Z
        ax.plot_surface(X2, Y2, Z2, color='orange', alpha=0.6, edgecolor='none')

    def plotar_paraboloide_hiperbolico(self, ax, params):
        # Implementação do parabolóide hiperbólico
        a = params.get('a', 1.0)
        b = params.get('b', 1.0)
        h = params.get('h', 0.0)
        k = params.get('k', 0.0)
        l = params.get('l', 0.0)
        rotacao_x = params.get('rotacao_x', 0.0)
        rotacao_y = params.get('rotacao_y', 0.0)
        rotacao_z = params.get('rotacao_z', 0.0)

        x = np.linspace(-5, 5, 100)
        y = np.linspace(-5, 5, 100)
        X, Y = np.meshgrid(x, y)
        Z = (X**2 / a**2) - (Y**2 / b**2) + l

        # Aplicar rotações (opcional)
        if rotacao_x != 0 or rotacao_y != 0 or rotacao_z != 0:
            # Rotações em ordem Z, Y, X
            rot_z = np.radians(rotacao_z)
            rot_y = np.radians(rotacao_y)
            rot_x = np.radians(rotacao_x)

            # Matriz de rotação Z
            Rz = np.array([
                [np.cos(rot_z), -np.sin(rot_z), 0],
                [np.sin(rot_z),  np.cos(rot_z), 0],
                [0,             0,             1]
            ])

            # Matriz de rotação Y
            Ry = np.array([
                [ np.cos(rot_y), 0, np.sin(rot_y)],
                [0,              1, 0],
                [-np.sin(rot_y), 0, np.cos(rot_y)]
            ])

            # Matriz de rotação X
            Rx = np.array([
                [1, 0,              0],
                [0, np.cos(rot_x), -np.sin(rot_x)],
                [0, np.sin(rot_x),  np.cos(rot_x)]
            ])

            # Aplicar rotações
            for i in range(X.shape[0]):
                for j in range(X.shape[1]):
                    vec = np.array([X[i, j], Y[i, j], Z[i, j]])
                    vec_rot = Rz @ Ry @ Rx @ vec
                    X[i, j], Y[i, j], Z[i, j] = vec_rot

        ax.plot_surface(X + h, Y + k, Z, color='purple', alpha=0.6, edgecolor='none', label='Parabolóide Hiperbólico')

    def plotar_paraboloide_eliptico(self, ax, params):
        # Implementação do parabolóide elíptico
        a = params.get('a', 1.0)
        b = params.get('b', 1.0)
        h = params.get('h', 0.0)
        k = params.get('k', 0.0)
        l = params.get('l', 0.0)
        rotacao_x = params.get('rotacao_x', 0.0)
        rotacao_y = params.get('rotacao_y', 0.0)
        rotacao_z = params.get('rotacao_z', 0.0)

        x = np.linspace(-5, 5, 100)
        y = np.linspace(-5, 5, 100)
        X, Y = np.meshgrid(x, y)
        Z = (X**2 / a**2) + (Y**2 / b**2) + l

        # Aplicar rotações (opcional)
        if rotacao_x != 0 or rotacao_y != 0 or rotacao_z != 0:
            # Rotações em ordem Z, Y, X
            rot_z = np.radians(rotacao_z)
            rot_y = np.radians(rotacao_y)
            rot_x = np.radians(rotacao_x)

            # Matriz de rotação Z
            Rz = np.array([
                [np.cos(rot_z), -np.sin(rot_z), 0],
                [np.sin(rot_z),  np.cos(rot_z), 0],
                [0,             0,             1]
            ])

            # Matriz de rotação Y
            Ry = np.array([
                [ np.cos(rot_y), 0, np.sin(rot_y)],
                [0,              1, 0],
                [-np.sin(rot_y), 0, np.cos(rot_y)]
            ])

            # Matriz de rotação X
            Rx = np.array([
                [1, 0,              0],
                [0, np.cos(rot_x), -np.sin(rot_x)],
                [0, np.sin(rot_x),  np.cos(rot_x)]
            ])

            # Aplicar rotações
            for i in range(X.shape[0]):
                for j in range(X.shape[1]):
                    vec = np.array([X[i, j], Y[i, j], Z[i, j]])
                    vec_rot = Rz @ Ry @ Rx @ vec
                    X[i, j], Y[i, j], Z[i, j] = vec_rot

        ax.plot_surface(X + h, Y + k, Z, color='pink', alpha=0.6, edgecolor='none', label='Parabolóide Elíptico')

    def plotar_superficie_conica(self, ax, params):
        # Implementação de uma superfície cônica
        a = params.get('a', 1.0)
        b = params.get('b', 1.0)
        c = params.get('c', 1.0)
        h = params.get('h', 0.0)
        k = params.get('k', 0.0)
        l = params.get('l', 0.0)
        rotacao_x = params.get('rotacao_x', 0.0)
        rotacao_y = params.get('rotacao_y', 0.0)
        rotacao_z = params.get('rotacao_z', 0.0)

        z = np.linspace(-5, 5, 100)
        r = np.linspace(0, 5, 100)
        R, Z = np.meshgrid(r, z)
        X = R * np.cos(Z) + h
        Y = R * np.sin(Z) + k
        Z = Z * c + l

        # Aplicar rotações (opcional)
        if rotacao_x != 0 or rotacao_y != 0 or rotacao_z != 0:
            # Rotações em ordem Z, Y, X
            rot_z = np.radians(rotacao_z)
            rot_y = np.radians(rotacao_y)
            rot_x = np.radians(rotacao_x)

            # Matriz de rotação Z
            Rz = np.array([
                [np.cos(rot_z), -np.sin(rot_z), 0],
                [np.sin(rot_z),  np.cos(rot_z), 0],
                [0,             0,             1]
            ])

            # Matriz de rotação Y
            Ry = np.array([
                [ np.cos(rot_y), 0, np.sin(rot_y)],
                [0,              1, 0],
                [-np.sin(rot_y), 0, np.cos(rot_y)]
            ])

            # Matriz de rotação X
            Rx = np.array([
                [1, 0,              0],
                [0, np.cos(rot_x), -np.sin(rot_x)],
                [0, np.sin(rot_x),  np.cos(rot_x)]
            ])

            # Aplicar rotações
            for i in range(X.shape[0]):
                for j in range(X.shape[1]):
                    vec = np.array([X[i, j], Y[i, j], Z[i, j]])
                    vec_rot = Rz @ Ry @ Rx @ vec
                    X[i, j], Y[i, j], Z[i, j] = vec_rot

        ax.plot_surface(X, Y, Z, color='brown', alpha=0.6, edgecolor='none', label='Superfície Cônica')

    def plotar_superficie_cilindrica(self, ax, params):
        # Implementação de uma superfície cilíndrica (cilindro)
        a = params.get('a', 1.0)  # Raio do cilindro
        h = params.get('h', 0.0)
        k = params.get('k', 0.0)
        l = params.get('l', 0.0)
        rotacao_x = params.get('rotacao_x', 0.0)
        rotacao_y = params.get('rotacao_y', 0.0)
        rotacao_z = params.get('rotacao_z', 0.0)

        z = np.linspace(-5, 5, 100)
        theta = np.linspace(0, 2 * np.pi, 100)
        Theta, Z = np.meshgrid(theta, z)
        X = a * np.cos(Theta) + h
        Y = a * np.sin(Theta) + k

        # Aplicar rotações (opcional)
        if rotacao_x != 0 or rotacao_y != 0 or rotacao_z != 0:
            # Rotações em ordem Z, Y, X
            rot_z = np.radians(rotacao_z)
            rot_y = np.radians(rotacao_y)
            rot_x = np.radians(rotacao_x)

            # Matriz de rotação Z
            Rz = np.array([
                [np.cos(rot_z), -np.sin(rot_z), 0],
                [np.sin(rot_z),  np.cos(rot_z), 0],
                [0,             0,             1]
            ])

            # Matriz de rotação Y
            Ry = np.array([
                [ np.cos(rot_y), 0, np.sin(rot_y)],
                [0,              1, 0],
                [-np.sin(rot_y), 0, np.cos(rot_y)]
            ])

            # Matriz de rotação X
            Rx = np.array([
                [1, 0,              0],
                [0, np.cos(rot_x), -np.sin(rot_x)],
                [0, np.sin(rot_x),  np.cos(rot_x)]
            ])

            # Aplicar rotações
            for i in range(X.shape[0]):
                for j in range(X.shape[1]):
                    vec = np.array([X[i, j], Y[i, j], Z[i, j]])
                    vec_rot = Rz @ Ry @ Rx @ vec
                    X[i, j], Y[i, j], Z[i, j] = vec_rot

        ax.plot_surface(X, Y, Z + l, color='grey', alpha=0.6, edgecolor='none', label='Superfície Cilíndrica')

    def plotar_superficie_esferica(self, ax, params):
        # Implementação de uma superfície esférica
        a = params.get('a', 1.0)  # Raio da esfera
        h = params.get('h', 0.0)
        k = params.get('k', 0.0)
        l = params.get('l', 0.0)
        rotacao_x = params.get('rotacao_x', 0.0)
        rotacao_y = params.get('rotacao_y', 0.0)
        rotacao_z = params.get('rotacao_z', 0.0)

        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        U, V = np.meshgrid(u, v)
        X = a * np.cos(U) * np.sin(V) + h
        Y = a * np.sin(U) * np.sin(V) + k
        Z = a * np.cos(V) + l

        # Aplicar rotações (opcional)
        if rotacao_x != 0 or rotacao_y != 0 or rotacao_z != 0:
            # Rotações em ordem Z, Y, X
            rot_z = np.radians(rotacao_z)
            rot_y = np.radians(rotacao_y)
            rot_x = np.radians(rotacao_x)

            # Matriz de rotação Z
            Rz = np.array([
                [np.cos(rot_z), -np.sin(rot_z), 0],
                [np.sin(rot_z),  np.cos(rot_z), 0],
                [0,             0,             1]
            ])

            # Matriz de rotação Y
            Ry = np.array([
                [ np.cos(rot_y), 0, np.sin(rot_y)],
                [0,              1, 0],
                [-np.sin(rot_y), 0, np.cos(rot_y)]
            ])

            # Matriz de rotação X
            Rx = np.array([
                [1, 0,              0],
                [0, np.cos(rot_x), -np.sin(rot_x)],
                [0, np.sin(rot_x),  np.cos(rot_x)]
            ])

            # Aplicar rotações
            for i in range(X.shape[0]):
                for j in range(X.shape[1]):
                    vec = np.array([X[i, j], Y[i, j], Z[i, j]])
                    vec_rot = Rz @ Ry @ Rx @ vec
                    X[i, j], Y[i, j], Z[i, j] = vec_rot

        ax.plot_surface(X, Y, Z, color='cyan', alpha=0.6, edgecolor='none', label='Superfície Esférica')

    # Funções para interações do mouse na plotagem 3D
    def on_mouse_press(self, event):
        self.arrastando = True
        self.ultimo_x = event.xdata if event.xdata else 0
        self.ultimo_y = event.ydata if event.ydata else 0

    def on_mouse_release(self, event):
        self.arrastando = False

    def on_mouse_move(self, event):
        if not self.arrastando:
            return

        if event.xdata and event.ydata:
            dx = event.xdata - self.ultimo_x
            dy = event.ydata - self.ultimo_y

            if event.button == 1:  # Rotação
                self.rotacao[0] += dy * 2
                self.rotacao[1] += dx * 2
            elif event.button == 3:  # Translação
                self.translacao[0] += dx
                self.translacao[1] += dy

            self.ultimo_x = event.xdata
            self.ultimo_y = event.ydata
            self.atualizar_plotagem()

    def on_scroll(self, event):
        # Zoom
        if event.button == 'up':
            self.translacao[2] += 0.1
        else:
            self.translacao[2] -= 0.1
        self.atualizar_plotagem()

    def configurar_pagina_mapa_conceitual(self):
        self.frame_mapa_conceitual = ttk.Frame(self.pagina_mapa_conceitual)
        self.frame_mapa_conceitual.pack(fill=tk.BOTH, expand=True)

        # Barra de ferramentas do mapa conceitual
        self.toolbar_mapa = ttk.Frame(self.frame_mapa_conceitual)
        self.toolbar_mapa.pack(fill=tk.X)

        ttk.Button(
            self.toolbar_mapa,
            text="Novo Conceito",
            command=self.adicionar_conceito
        ).pack(side=tk.LEFT, padx=5)
        ttk.Button(
            self.toolbar_mapa,
            text="Nova Relação",
            command=self.adicionar_relacao
        ).pack(side=tk.LEFT, padx=5)
        ttk.Button(
            self.toolbar_mapa,
            text="Salvar Mapa",
            command=self.salvar_mapa_conceitual
        ).pack(side=tk.LEFT, padx=5)
        ttk.Button(
            self.toolbar_mapa,
            text="Carregar Mapa",
            command=self.carregar_mapa_conceitual
        ).pack(side=tk.LEFT, padx=5)

        # Canvas para o mapa conceitual
        self.canvas_mapa = tk.Canvas(self.frame_mapa_conceitual, bg='white')
        self.canvas_mapa.pack(fill=tk.BOTH, expand=True)

        # Dicionários para armazenar conceitos e relações
        self.conceitos = {}
        self.relacoes = {}

        # Adicionar funcionalidade de arrastar
        self.adicionar_funcionalidade_arrastar()

    def adicionar_conceito(self):
        dialogo = tk.Toplevel(self.raiz)
        dialogo.title("Novo Conceito")

        ttk.Label(dialogo, text="Nome do Conceito:").pack(pady=5)
        nome_var = tk.StringVar()
        ttk.Entry(dialogo, textvariable=nome_var).pack(pady=5)

        ttk.Label(dialogo, text="Referências (separadas por vírgula):").pack(pady=5)
        referencias_var = tk.StringVar()
        ttk.Entry(dialogo, textvariable=referencias_var).pack(pady=5)

        def criar():
            nome = nome_var.get().strip()
            referencias = referencias_var.get().strip()
            if nome:
                x = len(self.conceitos) * 100 + 50
                y = 100
                id_conceito = self.canvas_mapa.create_oval(
                    x - 30, y - 20, x + 30, y + 20,
                    fill='lightblue'
                )
                id_texto = self.canvas_mapa.create_text(x, y, text=nome)
                id_ref = self.canvas_mapa.create_text(
                    x, y + 25,
                    text=f"Ref: {referencias}",
                    font=('Arial', 8),
                    fill='gray'
                )
                self.conceitos[nome] = {
                    'x': x,
                    'y': y,
                    'id': id_conceito,
                    'texto_id': id_texto,
                    'ref_id': id_ref,
                    'referencias': referencias  # Adiciona as referências
                }
            dialogo.destroy()

        ttk.Button(dialogo, text="Criar", command=criar).pack(pady=5)

    def adicionar_relacao(self):
        if len(self.conceitos) < 2:
            messagebox.showwarning("Aviso", "Precisa de pelo menos 2 conceitos!")
            return

        dialogo = tk.Toplevel(self.raiz)
        dialogo.title("Nova Relação")

        ttk.Label(dialogo, text="De:").pack(pady=5)
        de_var = tk.StringVar()
        ttk.Combobox(
            dialogo,
            textvariable=de_var,
            values=list(self.conceitos.keys()),
            state="readonly"
        ).pack(pady=5)

        ttk.Label(dialogo, text="Para:").pack(pady=5)
        para_var = tk.StringVar()
        ttk.Combobox(
            dialogo,
            textvariable=para_var,
            values=list(self.conceitos.keys()),
            state="readonly"
        ).pack(pady=5)

        ttk.Label(dialogo, text="Descrição:").pack(pady=5)
        desc_var = tk.StringVar()
        ttk.Entry(dialogo, textvariable=desc_var).pack(pady=5)

        def criar():
            de = de_var.get()
            para = para_var.get()
            descricao = desc_var.get().strip()

            if de and para and descricao:
                x1 = self.conceitos[de]['x']
                y1 = self.conceitos[de]['y']
                x2 = self.conceitos[para]['x']
                y2 = self.conceitos[para]['y']

                id_linha = self.canvas_mapa.create_line(
                    x1, y1, x2, y2, arrow=tk.LAST
                )
                id_texto = self.canvas_mapa.create_text(
                    (x1 + x2) / 2,
                    (y1 + y2) / 2,
                    text=descricao,
                    font=('Arial', 10),
                    fill='black'
                )

                chave_relacao = f"{de}-{para}-{descricao}"
                self.relacoes[chave_relacao] = {
                    'de': de,
                    'para': para,
                    'descricao': descricao,
                    'linha_id': id_linha,
                    'texto_id': id_texto
                }
            dialogo.destroy()

        ttk.Button(dialogo, text="Criar", command=criar).pack(pady=5)

    def salvar_mapa_conceitual(self):
        dados = {
            'conceitos': self.conceitos,
            'relacoes': self.relacoes
        }
        try:
            with open('mapa_conceitual.json', 'w') as f:
                json.dump(dados, f, indent=4)
            messagebox.showinfo("Sucesso", "Mapa conceitual salvo com sucesso!")
        except Exception as e:
            messagebox.showerror("Erro", f"Erro ao salvar mapa: {str(e)}")

    def carregar_mapa_conceitual(self):
        try:
            with open('mapa_conceitual.json', 'r') as f:
                dados = json.load(f)

            # Limpar canvas
            self.canvas_mapa.delete('all')
            self.conceitos.clear()
            self.relacoes.clear()

            # Recriar conceitos
            self.conceitos = dados.get('conceitos', {})
            for nome, dados_conceito in self.conceitos.items():
                x = dados_conceito['x']
                y = dados_conceito['y']
                referencias = dados_conceito.get('referencias', '')
                id_conceito = self.canvas_mapa.create_oval(
                    x - 30, y - 20, x + 30, y + 20,
                    fill='lightblue'
                )
                id_texto = self.canvas_mapa.create_text(x, y, text=nome)
                id_ref = self.canvas_mapa.create_text(
                    x, y + 25,
                    text=f"Ref: {referencias}",
                    font=('Arial', 8),
                    fill='gray'
                )
                self.conceitos[nome]['id'] = id_conceito
                self.conceitos[nome]['texto_id'] = id_texto
                self.conceitos[nome]['ref_id'] = id_ref

            # Recriar relações
            self.relacoes = dados.get('relacoes', {})
            for chave, dados_relacao in self.relacoes.items():
                de = dados_relacao['de']
                para = dados_relacao['para']
                descricao = dados_relacao['descricao']
                x1 = self.conceitos[de]['x']
                y1 = self.conceitos[de]['y']
                x2 = self.conceitos[para]['x']
                y2 = self.conceitos[para]['y']

                id_linha = self.canvas_mapa.create_line(
                    x1, y1, x2, y2, arrow=tk.LAST
                )
                id_texto = self.canvas_mapa.create_text(
                    (x1 + x2) / 2,
                    (y1 + y2) / 2,
                    text=descricao,
                    font=('Arial', 10),
                    fill='black'
                )

                self.relacoes[chave]['linha_id'] = id_linha
                self.relacoes[chave]['texto_id'] = id_texto

            messagebox.showinfo("Sucesso", "Mapa conceitual carregado com sucesso!")
        except FileNotFoundError:
            messagebox.showerror("Erro", "Arquivo de mapa conceitual não encontrado!")
        except Exception as e:
            messagebox.showerror("Erro", f"Erro ao carregar mapa: {str(e)}")

    def adicionar_funcionalidade_arrastar(self):
        self.canvas_mapa.bind('<Button-1>', self.iniciar_arrasto)
        self.canvas_mapa.bind('<B1-Motion>', self.arrastando_conceito)
        self.canvas_mapa.bind('<ButtonRelease-1>', self.parar_arrasto)
        self.item_arrastando = None
        self.deslocamento_x = 0
        self.deslocamento_y = 0

    def iniciar_arrasto(self, event):
        # Encontrar o item mais próximo ao clique
        item = self.canvas_mapa.find_closest(event.x, event.y)
        if item:
            # Verificar se o item é um conceito (oval) baseado nas coordenadas
            coords = self.canvas_mapa.coords(item)
            if len(coords) == 4:  # Oval tem 4 coordenadas
                self.item_arrastando = item[0]
                self.deslocamento_x = event.x
                self.deslocamento_y = event.y

    def arrastando_conceito(self, event):
        if self.item_arrastando:
            dx = event.x - self.deslocamento_x
            dy = event.y - self.deslocamento_y
            self.canvas_mapa.move(self.item_arrastando, dx, dy)
            self.deslocamento_x = event.x
            self.deslocamento_y = event.y

            # Atualizar posição no dicionário
            for nome, dados_conceito in self.conceitos.items():
                if dados_conceito['id'] == self.item_arrastando:
                    self.conceitos[nome]['x'] += dx
                    self.conceitos[nome]['y'] += dy
                    # Mover o texto junto
                    self.canvas_mapa.move(dados_conceito['texto_id'], dx, dy)
                    # Mover as referências
                    self.canvas_mapa.move(dados_conceito['ref_id'], dx, dy)
                    break

            # Atualizar relações conectadas
            self.atualizar_relacoes()

    def parar_arrasto(self, event):
        self.item_arrastando = None

    def atualizar_relacoes(self):
        """Atualiza as posições das linhas de relação quando os conceitos são movidos"""
        for relacao in self.relacoes.values():
            de = relacao['de']
            para = relacao['para']
            x1 = self.conceitos[de]['x']
            y1 = self.conceitos[de]['y']
            x2 = self.conceitos[para]['x']
            y2 = self.conceitos[para]['y']

            self.canvas_mapa.coords(relacao['linha_id'], x1, y1, x2, y2)

            # Atualizar posição do texto da relação
            self.canvas_mapa.coords(
                relacao['texto_id'],
                (x1 + x2) / 2,
                (y1 + y2) / 2
            )

if __name__ == "__main__":
    # Configurar estilo visual
    if TEMA_DISPONIVEL:
        raiz = ThemedTk(theme="arc")  # Tema moderno
    else:
        raiz = tk.Tk()

    # Configurar estilo padrão
    estilo = ttk.Style()
    estilo.configure("TButton", padding=6)
    estilo.configure("TLabel", padding=3)
    estilo.configure("TFrame", padding=5)
    estilo.configure("TEntry", padding=5)
    estilo.configure("TCombobox", padding=5)

    # Inicializar aplicação
    try:
        app = AplicativoVisualizadorGeometriaAvancada(raiz)

        # Configurar ícone (opcional)
        try:
            raiz.iconbitmap('icon.ico')  # Substitua pelo caminho do seu ícone
        except:
            pass

        # Iniciar loop principal
        raiz.mainloop()
    except Exception as e:
        messagebox.showerror("Erro", f"Erro ao iniciar aplicação: {str(e)}")
