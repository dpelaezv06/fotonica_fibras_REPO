import numpy as np
from scipy.optimize import minimize


def modos_TE(n_core, n_cleavy, espesor, angulo, modo, n_substract = None, longitud_onda = 1):
    ''' funcion que realiza el calculo de las fases y la ecuacion trascendednte para modos TE
    
    Entradas:
    n_core (float) == indice de refraccion del core del guia de onda
    n_cleavy (float) == indice de refraccion del recubrimiento del guia de onda
    espesor (float) == ancho del core del guia de onda (en micras)
    angulo (float) == angulo con el cual el rayo hace el zig-zag en el guia de onda
    modo (int) == numero del modo que se esta intentando evaluar en el guia de onda
    n_substract (float) == indice de refraccion del sustrato del guia de onda, por defecto es el mismo valor que se pone en el reccubrimiento
    longitud_onda (float) == longitud de onda de la iluminacion incidente en micras (por defecto es una micra)
    
    Retorna: Valor de la ecuacion trascendente usando trazado de rayos (esta funcion esta hecha para ser usada con minimize, los 
    angulos de entrada son solamente puntos de partida desde los cuales inicia la minimizacion)'''


    if n_substract is None: #en caso de que no se especifique el indice de refraccion del substract
        n_substract = n_cleavy #se toma como un guia de onda simetrico

    ''' condiciones de iluminacion '''
    numero_onda = 2 * np.pi / longitud_onda #numero de onda en el vacio k_0
    angulo = angulo * np.pi / 180 #un angulo inicial para realizar la propagacion

    ''' parametros de fase '''    
    fase_propagacion = n_core * numero_onda * espesor * np.cos(angulo) #fase acumulada por propagacion
    fase_reflexion = np.arctan(n_cleavy / (n_core * np.cos(angulo)) * np.sqrt((n_cleavy**2) / (n_core**2) * (np.sin(angulo))**2 - 1)) #fase por reflexion en las superficies
    fase_acumulada = 2 * (fase_propagacion - 2 * fase_reflexion) #fase acumulada en un zig-zag
    ecuacion_trascendente = fase_acumulada - modo * 2 * np.pi #ecuacion trascendente que habria que minimizar
    ecuacion_valorAbsoluto = np.abs(ecuacion_trascendente) #sacamos el valor absoluto parar evitar los valores negativos, y hacer que la ecuacion trascendente tienda a cero cuando se minimice
    return ecuacion_valorAbsoluto #retornamos el valor de la ecuacion tracendente, hay que hacer que sea cero para calcular el angulo de los modos propagantes
    

def modos_TM(n_core, n_cleavy, espesor, angulo, modo, n_substract = None, longitud_onda = 1):
    ''' funcion que realiza el calculo de las fases y la ecuacion trascendednte para modos TM
    
    Entradas:
    n_core (float) == indice de refraccion del core del guia de onda
    n_cleavy (float) == indice de refraccion del recubrimiento del guia de onda
    espesor (float) == ancho del core del guia de onda (en micras)
    angulo (float) == angulo con el cual el rayo hace el zig-zag en el guia de onda
    modo (int) == numero del modo que se esta intentando evaluar en el guia de onda
    n_substract (float) == indice de refraccion del sustrato del guia de onda, por defecto es el mismo valor que se pone en el reccubrimiento
    longitud_onda (float) == longitud de onda de la iluminacion incidente en micras (por defecto es una micra)
    
    Retorna: Valor de la ecuacion trascendente usando trazado de rayos (esta funcion esta hecha para ser usada con minimize, los 
    angulos de entrada son solamente puntos de partida desde los cuales inicia la minimizacion)'''


    if n_substract is None: #en caso de que no se especifique el indice de refraccion del substract
        n_substract = n_cleavy #se toma como un guia de onda simetrico

    ''' condiciones de iluminacion '''
    numero_onda = 2 * np.pi / longitud_onda #numero de onda en el vacio k_0
    angulo = angulo * np.pi / 180 #un angulo inicial para realizar la propagacion

    ''' parametros de fase '''    
    fase_propagacion = n_core * numero_onda * espesor * np.cos(angulo) #fase acumulada por propagacion
    fase_reflexion = np.arctan(n_core / (n_cleavy * np.cos(angulo)) * np.sqrt((n_cleavy**2) / (n_core**2) * (np.sin(angulo))**2 - 1)) #fase por reflexion en las superficies
    fase_acumulada = 2 * (fase_propagacion - 2 * fase_reflexion) #fase acumulada en un zig-zag
    ecuacion_trascendente = fase_acumulada - modo * 2 * np.pi #ecuacion trascendente que habria que minimizar
    ecuacion_valorAbsoluto = np.abs(ecuacion_trascendente) #sacamos el valor absoluto parar evitar los valores negativos, y hacer que la ecuacion trascendente tienda a cero cuando se minimice
    return ecuacion_valorAbsoluto #retornamos el valor de la ecuacion tracendente, hay que hacer que sea cero para calcular el angulo de los modos propagantes

def optimizar_TE(n_core, n_cleavy, espesor, modo, n_substract=None, longitud_onda=1):
    ''' funcion que calcula el angulo optimo para que la ecuacion trascendente de modos TE sea cero
    
    Entradas:
    n_core (float) == indice de refraccion del nucleo del guia de onda
    n_cleavy (float) == indice de refraccion del recubrimiento del guia de onda
    espesor (float) == ancho del nucleo del guia de onda (en micras)
    modo (int) == numero del modo que se esta intentando evaluar en el guia de onda
    n_substract (float, opcional) == indice de refraccion del sustrato (por defecto es el mismo valor que el recubrimiento)
    longitud_onda (float, opcional) == longitud de onda de la iluminacion incidente en micras (por defecto es una micra)
    
    Retorna:
    angulo_optimo (float) == angulo que hace que la ecuacion trascendente sea aproximadamente cero (en grados)
    valor_min (float) == valor minimo obtenido en la minimizacion (debe ser cercano a cero si la optimizacion es exitosa) '''

    # Valor inicial del angulo para comenzar la minimizacion (en grados)
    angulo_inicial = 45

    # Ejecucion de la minimizacion utilizando la funcion modos_TE
    resultado = minimize(
        modos_TE,                   # Funcion a minimizar (modos_TE)
        angulo_inicial,             # Punto de partida para el angulo (45 grados)
        args=(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda, True),  # Parametros adicionales de la funcion
        bounds=[(0, 90)]            # Restriccion del angulo entre 0 y 90 grados
    )
    
    # Retornar el angulo optimo encontrado y el valor minimo alcanzado
    return resultado.x[0], resultado.fun

def optimizar_TM(n_core, n_cleavy, espesor, modo, n_substract=None, longitud_onda=1):
    ''' funcion que calcula el angulo optimo para que la ecuacion trascendente de modos TM sea cero
    
    Entradas:
    n_core (float) == indice de refraccion del nucleo del guia de onda
    n_cleavy (float) == indice de refraccion del recubrimiento del guia de onda
    espesor (float) == ancho del nucleo del guia de onda (en micras)
    modo (int) == numero del modo que se esta intentando evaluar en el guia de onda
    n_substract (float, opcional) == indice de refraccion del sustrato (por defecto es el mismo valor que el recubrimiento)
    longitud_onda (float, opcional) == longitud de onda de la iluminacion incidente en micras (por defecto es una micra)
    
    Retorna:
    angulo_optimo (float) == angulo que hace que la ecuacion trascendente sea aproximadamente cero (en grados)
    valor_min (float) == valor minimo obtenido en la minimizacion (debe ser cercano a cero si la optimizacion es exitosa) '''
    
    angulo_critico = np.arcsin(n_cleavy/n_core) #se calcula el angulo critico del guia de onda
    angulo_inicial = 1.001*angulo_critico #se inicia con un valor de angulo critico mas el 0.1% del valor del angulo critico
    angulo_inicial = angulo_inicial * (np.pi/180) #se pasa el angulo inicial a grados

    ''' ejecucion de la minimizacion utilizando la funcion modos_TM
    funcion a minimizar (modos_TM)
    punto de partida para el angulo el angulo critico + el 0.01% del angulo critico
    parametros adicionales de la funcion
    restriccion del angulo entre el angulo critico y 90 grados '''
    resultado = minimize(modos_TM, angulo_inicial, args=(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda, True), bounds=[(angulo_critico, 90)])
    
    # Retornar el angulo optimo encontrado y el valor minimo alcanzado
    return resultado.x[0], resultado.fun

