import numpy as np
import modos_opticaRayos as ray

''' parametros del guia de onda '''
n_core = 1.5 #indice de refraccion del nucleo
n_cleavy = 1 #indice de refraccion del recubrimiento
n_substract = 1 #indice de refraccion del sustrato
espesor = 1 #espesor del guia de onda
longitud_onda = 1 #longitud de onda usada en la fibra optica

''' calculo de los valores solucion de la ecuacion trascendente'''
for modo in range(0, 4): #se recorren los valores de modos desde el 0 hasta el 3
    angulo_optimo, valor_min = ray.optimizar_TE(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda) #calculo del angulo para el modo TM
    ''' printeamos por consola el valor del angulo y de la ecuacion trascendednte en ese angulo. Si el valor de la ecuacion no es cero, entonces ese modo no es 
    permitido en el guia de onda '''
    print(f"Modo TE {modo}: Angulo optimo = {angulo_optimo*(180/np.pi):.6f}°, Valor minimo = {valor_min:.3e}") #se arroja el print de los datos


    angulo_optimo, valor_min = ray.optimizar_TM(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda) #calculo del angulo para el modo TE
    ''' printeamos por consola el valor del angulo y de la ecuacion trascendednte en ese angulo. Si el valor de la ecuacion no es cero, entonces ese modo no es 
    permitido en el guia de onda '''
    print(f"Modo TM {modo}: Angulo optimo = {angulo_optimo*(180/np.pi):.6f}°, Valor minimo = {valor_min:.3e}") #se arroja el print de los datos

