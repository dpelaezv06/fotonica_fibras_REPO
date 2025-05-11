import numpy as np
from scipy.optimize import minimize

def modos_TERayos(angulo, n_core, n_cleavy, espesor, modo, n_substract, longitud_onda):
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

    ''' condiciones de iluminacion '''
    numero_onda = 2 * np.pi / longitud_onda #numero de onda en el vacio k_0

    ''' parametros de fase '''    
    fase_propagacion = n_core * numero_onda * espesor * np.cos(angulo) #fase acumulada por propagacion
    fase_reflexion = np.arctan(n_cleavy / (n_core * np.cos(angulo)) * np.sqrt((n_core**2) * (np.sin(angulo))**2 / (n_cleavy**2) - 1)) #fase por reflexion en las superficies
    fase_acumulada = 2 * (fase_propagacion - 2 * fase_reflexion) #fase acumulada en un zig-zag
    ecuacion_trascendente = fase_acumulada - modo * 2 * np.pi #ecuacion trascendente que habria que minimizar
    ecuacion_valorAbsoluto = np.abs(ecuacion_trascendente) #sacamos el valor absoluto parar evitar los valores negativos, y hacer que la ecuacion trascendente tienda a cero cuando se minimice
    return ecuacion_valorAbsoluto #retornamos el valor de la ecuacion tracendente, hay que hacer que sea cero para calcular el angulo de los modos propagantes
    

def modos_TMRayos(angulo, n_core, n_cleavy, espesor, modo, n_substract, longitud_onda):
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
    ''' condiciones de iluminacion '''
    numero_onda = 2 * np.pi / longitud_onda #numero de onda en el vacio k_0

    ''' parametros de fase '''    
    fase_propagacion = n_core * numero_onda * espesor * np.cos(angulo) #fase acumulada por propagacion
    fase_reflexion = np.arctan(n_core / (n_cleavy * np.cos(angulo)) * np.sqrt(((n_core**2) / (n_cleavy**2) * (np.sin(angulo))**2) - 1)) #fase por reflexion en las superficies
    fase_acumulada = 2 * (fase_propagacion - 2 * fase_reflexion) #fase acumulada en un zig-zag
    ecuacion_trascendente = fase_acumulada - modo * 2 * np.pi #ecuacion trascendente que habria que minimizar
    ecuacion_valorAbsoluto = np.abs(ecuacion_trascendente) #sacamos el valor absoluto parar evitar los valores negativos, y hacer que la ecuacion trascendente tienda a cero cuando se minimice
    return ecuacion_valorAbsoluto #retornamos el valor de la ecuacion tracendente, hay que hacer que sea cero para calcular el angulo de los modos propagantes

def optimizar_TERayos(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda):
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

    angulo_critico = np.arcsin(n_cleavy/n_core) #se calcula el angulo critico del guia de onda
    angulo_inicial = 1.0001*angulo_critico #se inicia con un valor de angulo critico mas el 0.01% del valor del angulo critico

    ''' ejecucion de la minimizacion utilizando la funcion modos_TE, el punto de partida para el angulo es el angulo critico + el 0.01% del angulo critico,
    se aplica una restriccion del angulo entre el angulo critico y 90 grados '''
    resultado = minimize(modos_TERayos, angulo_inicial, args=(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda), bounds=[(angulo_critico, np.pi/2)])
    
    # Retornar el angulo optimo encontrado y el valor minimo alcanzado
    return resultado.x[0], resultado.fun

def optimizar_TMRayos(n_core, n_cleavy, espesor, modo, n_substract=None, longitud_onda=1):
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
    angulo_inicial = 1.0001*angulo_critico #se inicia con un valor de angulo critico mas el 0.01% del valor del angulo critico

    ''' ejecucion de la minimizacion utilizando la funcion modos_TM, el punto de partida para el angulo es el angulo critico + el 0.01% del angulo critico,
    se aplica una restriccion del angulo entre el angulo critico y 90 grados '''
    resultado = minimize(modos_TMRayos, angulo_inicial, args=(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda), bounds=[(angulo_critico, np.pi/2)])
    
    # Retornar el angulo optimo encontrado y el valor minimo alcanzado
    return resultado.x[0], resultado.fun

def n_efectivo(n_core, angulo):
    ''' funcion que devuelve los parametros de indice de refraccion efectivo de un modo, dependiendo del angulo de incidencia dentro del guia de onda y el 
    indice de refraccion del nucleo del guia de onda
    
    ENTRADAS:
    n_core (float) == indice de refraccion del nucleo del guia de onda
    angulo (float) == angulo (EN RADIANES) de incidencia del zig zag del rayo dentro del guia de onda, es relativo a cada modo de propagacion '''
    n_efectivo = n_core*np.sin(angulo)
    return n_efectivo
