import numpy as np
from scipy.optimize import minimize

def modos_TEOndasPares(angulo, modo, n_core, n_cleavy, espesor, n_substract, longitud_onda):
    ''' ecuacion de los modos pares para los modos TE calculados a partir de la teoria ondulatoria, retorna el valor absoluto de la ecuacion trascendente
    para ser minimizado
    ENTRADAS:
    angulo (float) == angulo del zigzag de la luz en el guia de onda
    modo (int) == modo del guia de onda que se va a examinar
    n_core (float) == indice de reefraccion del core del guia de onda
    n_cleavy (float) == indice de refrarccion del cleavy del guia de onda
    espesor (float) == espesor del guia de onda en micras
    n_substract (float) == indice de refraccion del substract del guia de onda
    longitud_onda (float) == longitud de onda de la iluminacion incidente en el guia de onda, en micras 
    
    RETORNA:
    Valor absoluto del resultado de la ecuacion trascendente despejado a cero, con el fin de minimizar '''

    numero_onda = 2 * np.pi / longitud_onda #numero de onda en el vacio de la iluminacion del guia de onda
    n_efectivo = n_core * np.sin(angulo) #calculo del indice de refraccion efectivo
    gamma = numero_onda * np.sqrt(n_efectivo**2 - n_cleavy**2) #calculo del factor gamma, el cual depende de las condiciones del recubrimiento
    kappa = numero_onda * np.sqrt(n_core**2 - n_efectivo**2) #calculo del factor kappa, el cual depende de las condiciones del core
    ecuacion_trascendente = kappa * espesor/2 - np.arctan(gamma/kappa) + modo*np.pi #ecuacion trascendente igualada a cero
    ecuacion_paresAbsoluto = np.abs(ecuacion_trascendente) #se calcula el valor absoluto para poder optimizar con minimizacion tendiendo a cero
    return ecuacion_paresAbsoluto #se retorna el valor absoluto del valor de la ecuacion trascendente

def modos_TEOndasImpares(angulo, modo, n_core, n_cleavy, espesor, n_substract, longitud_onda):
    ''' ecuacion de los modos impares para los modos TE calculados a partir de la teoria ondulatoria, retorna el valor absoluto de la ecuacion trascendente
    para ser minimizado
    ENTRADAS:
    angulo (float) == angulo del zigzag de la luz en el guia de onda
    modo (int) == modo del guia de onda que se va a examinar
    n_core (float) == indice de reefraccion del core del guia de onda
    n_cleavy (float) == indice de refrarccion del cleavy del guia de onda
    espesor (float) == espesor del guia de onda en micras
    n_substract (float) == indice de refraccion del substract del guia de onda
    longitud_onda (float) == longitud de onda de la iluminacion incidente en el guia de onda, en micras 
    
    RETORNA:
    Valor absoluto del resultado de la ecuacion trascendente despejado a cero, con el fin de minimizar '''

    numero_onda = 2 * np.pi / longitud_onda #numero de onda de la iluminacion del guia de onda
    n_efectivo = n_core * np.sin(angulo) #calculo del indice de refraccion efectivo
    gamma = numero_onda * np.sqrt(n_efectivo**2 - n_cleavy**2) #calculo del factor gamma, el cual depende de las condiciones del recubrimiento
    kappa = numero_onda * np.sqrt(n_core**2 - n_efectivo**2) #calculo del factor kappa, el cual depende de las condiciones del core
    ecuacion_trascendente = kappa * espesor/2 - np.arctan(-kappa/gamma) + modo*np.pi #calculo del valor de la ecuacion trascendente despejada a cero
    ecuacion_imparesAbsoluto = np.abs(ecuacion_trascendente) #se calcula el valor absoluto de la funcion despejada a cero para que minimize funcione a cero
    return ecuacion_imparesAbsoluto #se retorna el valor de la ecuacion trascendente despejado a cero en valor absoluto para minimizarlo

def optimizar_TEOndasPares(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda):
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
    resultado = minimize(modos_TEOndasPares, angulo_inicial, args=(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda), bounds=[(angulo_critico, np.pi/2)])
    
    # Retornar el angulo optimo encontrado y el valor minimo alcanzado
    return resultado.x[0], resultado.fun

def optimizar_TEOndasImPares(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda):
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
    resultado = minimize(modos_TEOndasImpares, angulo_inicial, args=(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda), bounds=[(angulo_critico, np.pi/2)])
    
    # Retornar el angulo optimo encontrado y el valor minimo alcanzado
    return resultado.x[0], resultado.fun

def modos_TMOndasPares(angulo, modo, n_core, n_cleavy, espesor, n_substract, longitud_onda):
    ''' ecuacion de los modos pares para los modos TM calculados a partir de la teoria ondulatoria, retorna el valor absoluto de la ecuacion trascendente
    para ser minimizado
    ENTRADAS:
    angulo (float) == angulo del zigzag de la luz en el guia de onda
    modo (int) == modo del guia de onda que se va a examinar
    n_core (float) == indice de reefraccion del core del guia de onda
    n_cleavy (float) == indice de refrarccion del cleavy del guia de onda
    espesor (float) == espesor del guia de onda en micras
    n_substract (float) == indice de refraccion del substract del guia de onda
    longitud_onda (float) == longitud de onda de la iluminacion incidente en el guia de onda, en micras 
    
    RETORNA:
    Valor absoluto del resultado de la ecuacion trascendente despejado a cero, con el fin de minimizar '''

    numero_onda = 2*np.pi / longitud_onda #calculo del numero de onda en el vacio k_0
    n_efectivo = n_core * np.sin(angulo) #calculo del indice de refraccion efectivo
    kappa = numero_onda * np.sqrt(n_core**2 - n_efectivo**2) #calculo del parametro kappa que esta relacionado con la onda en el core
    gamma = numero_onda * np.sqrt(n_efectivo**2 - n_cleavy**2) #calculo del parametro gamma que esta relacionado con la onda en el cleavy
    tangente = np.arctan((gamma*n_core**2)/kappa) #termino de la tangente inversa en la ecuacio trascendente
    ecuacion_trascendente = kappa * espesor/2 + modo*np.pi - tangente #ecuacion trascendente, se debe resolver, esta igualada a cero
    ecuacion_trascendenteValorAbsoluto = np.abs(ecuacion_trascendente) #se saca el valor absoluto para que el minimo valor que tome sea cero y se pueda minimizar
    return ecuacion_trascendenteValorAbsoluto #se retorna el valor absoluto de la ecuacion trascendente, para pooder minimzar esta funcion

def modos_TMOndasImpares(angulo, modo, n_core, n_cleavy, espesor, n_substract, longitud_onda):
    ''' ecuacion de los modos impares para los modos TM calculados a partir de la teoria ondulatoria, retorna el valor absoluto de la ecuacion trascendente
    para ser minimizado
    ENTRADAS:
    angulo (float) == angulo del zigzag de la luz en el guia de onda
    modo (int) == modo del guia de onda que se va a examinar
    n_core (float) == indice de reefraccion del core del guia de onda
    n_cleavy (float) == indice de refrarccion del cleavy del guia de onda
    espesor (float) == espesor del guia de onda en micras
    n_substract (float) == indice de refraccion del substract del guia de onda
    longitud_onda (float) == longitud de onda de la iluminacion incidente en el guia de onda, en micras 
    
    RETORNA:
    Valor absoluto del resultado de la ecuacion trascendente despejado a cero, con el fin de minimizar '''


    numero_onda = 2*np.pi/longitud_onda #calculo del numero de onda en el vacio k_0
    n_efectivo = n_core * np.sin(angulo) #calculo del indice de refraccion efectivo
    kappa = numero_onda * np.sqrt(n_core**2 - n_efectivo**2) #calculo del parametro kappa que esta relacionado con la onda en el core
    gamma = numero_onda * np.sqrt(n_efectivo**2 - n_cleavy**2) #calculo del parametro gamma que esta relacionado con la onda en el cleavy
    tangente = np.arctan((-kappa) / (gamma * n_core**2)) #termino de la tangente inversa en la ecuacio trascendente
    ecuacion_trascendente = kappa * espesor/2 + modo*np.pi - tangente #ecuacion trascendente, se debe resolver, esta igualada a cero
    ecuacion_trascendenteValorAbsoluto = np.abs(ecuacion_trascendente) #se saca el valor absoluto para que el minimo valor que tome sea cero y se pueda minimizar
    return ecuacion_trascendenteValorAbsoluto #se retorna el valor absoluto de la ecuacion trascendente, para pooder minimzar esta funcion

def optimizar_TMOndasPares(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda):
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
    resultado = minimize(modos_TMOndasPares, angulo_inicial, args=(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda), bounds=[(angulo_critico, np.pi/2)])
    
    # Retornar el angulo optimo encontrado y el valor minimo alcanzado
    return resultado.x[0], resultado.fun

def optimizar_TMOndasImPares(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda):
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
    resultado = minimize(modos_TMOndasImpares, angulo_inicial, args=(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda), bounds=[(angulo_critico, np.pi/2)])
    
    # Retornar el angulo optimo encontrado y el valor minimo alcanzado
    return resultado.x[0], resultado.fun