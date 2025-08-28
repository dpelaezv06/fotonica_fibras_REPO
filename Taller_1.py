import numpy as np
import modos_opticaRayos as ray
import modos_opticaOndulatoria as ond


def n_efectivo(n_core, angulo):
    ''' funcion que devuelve los parametros de indice de refraccion efectivo de un modo, dependiendo del angulo de incidencia dentro del guia de onda y el 
    indice de refraccion del nucleo del guia de onda
    
    ENTRADAS:
    n_core (float) == indice de refraccion del nucleo del guia de onda
    angulo (float) == angulo (EN RADIANES) de incidencia del zig zag del rayo dentro del guia de onda, es relativo a cada modo de propagacion '''
    n_efectivo = n_core*np.sin(angulo)
    return n_efectivo

''' parametros del guia de onda '''
n_core = 1.5 #indice de refraccion del nucleo
n_cleavy = 1 #indice de refraccion del recubrimiento
n_substract = 1 #indice de refraccion del sustrato
espesor = 1 #espesor del guia de onda
longitud_onda = 1 #longitud de onda usada en la fibra optica

''' calculo de los valores solucion de la ecuacion trascendente'''
print("OPTICA DE RAYOS \n")
for modo in range(0, 3): #se recorren los valores de modos desde el 0 hasta el 3
    angulo_optimo, valor_min = ray.optimizar_TERayos(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda) #calculo del angulo para el modo TM
    ''' printeamos por consola el valor del angulo y de la ecuacion trascendednte en ese angulo. Si el valor de la ecuacion no es cero, entonces ese modo no es 
    permitido en el guia de onda '''
    indice_efectivo = n_efectivo(n_core, angulo_optimo)
    print(f"Modo TE {modo}: Angulo optimo = {angulo_optimo:.6f}, n efectivo = {indice_efectivo:.6f}") #se arroja el print de los datos


    angulo_optimo, valor_min = ray.optimizar_TMRayos(n_core, n_cleavy, espesor, modo, n_substract, longitud_onda) #calculo del angulo para el modo TE
    indice_efectivo = n_efectivo(n_core, angulo_optimo)
    ''' printeamos por consola el valor del angulo y de la ecuacion trascendednte en ese angulo. Si el valor de la ecuacion no es cero, entonces ese modo no es 
    permitido en el guia de onda '''
    print(f"Modo TM {modo}: Angulo optimo = {angulo_optimo:.6f}, n efectivo = {indice_efectivo:.6f}") #se arroja el print de los datos

print("\n OPTICA ONDULATORIA \n")

for modo in range(1): #se recorren los valores de modos desde el 0 hasta el 3
    angulo_optimoPar, valor_minPar = ond.optimizar_TEOndasPares(modo, n_core, n_cleavy, espesor, n_substract, longitud_onda) #calculo del angulo para el modo TM
    indice_efectivoPar = n_efectivo(n_core, angulo_optimoPar)
    print(f"Modo TE Par {modo}: Angulo optimo = {angulo_optimoPar:.6f}, n efectivo = {indice_efectivoPar:.6f}") #se arroja el print de los datos
