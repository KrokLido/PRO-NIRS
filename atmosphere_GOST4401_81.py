"""
Черноморец Илья, 2022.
Данный модуль предоставляет вывод параметров атмосферы.
"""
import math

# Функция взятия высоты из файла *.csv.
def f_rho(altitude, fcsv=True):    #От 0м до 5 км - каждые 100 метров. После - каждые 500м. Предельная высота - 69.5 км.
    #Считаем весь файл в список.
    l_alt=[]
    l_rho=[]
    if(fcsv):
        f_out = open('atmosphere_GOST4401_81.csv', 'r')
        for line in f_out:
            line = line.replace('\n', '')
            line = line.replace(':', '')
            lwr=line.split(',')
            l_alt.append(float(lwr[0]))
            l_rho.append(float(lwr[1]))
    else:
        f_out = open('atmosphere_GOST4401_81_data_rho.txt', 'r')
        for line in f_out:
            line = line.replace('\n', '')
            lwr = line.split(',')
            l_alt.append(float(lwr[0]))
            l_rho.append(float(lwr[1]))
    
    if altitude <= 5000:
        ost = math.fmod(altitude, 100)
        if ost>=50:
            altitude_double = altitude + 100 - ost
            index = int(altitude_double/100)
            ret = float(l_rho[index])
            return (ret)
        else:
            altitude_double = altitude - ost
            index = int(altitude_double/100)
            ret = float(l_rho[index])
            return (ret)
    
    if(altitude >=5000 and altitude <= 51000):
        ost = math.fmod(altitude, 500)
        
        if ost < 250:
            altitude_double = altitude - math.fmod(altitude, 500)
            index = 51 + int(altitude_double/500) -3
            return(float(l_rho[index]))
        if ost >= 250:
            altitude_double = altitude - math.fmod(altitude, 500) + 500
            index = 51 + int(altitude_double/500) -3
            return(float(l_rho[index]))
    
    
    if(altitude >= 51000 and altitude < 69500):
        ost = math.fmod(altitude, 1000)
        if(ost < 500):
            altitude_double = altitude - math.fmod(altitude, 1000)
            index = 150 + int(altitude_double/1000) -51
            return(float(l_rho[index]))
        if( ost >= 500):
            altitude_double = altitude - math.fmod(altitude, 1000) + 1000
            index = 150 + int(altitude_double/1000) - 51
            return(float(l_rho[index]))



# Процедура получает словарь списков. Из них проводит выборку.
# Получает:				словарь списков. l_alt - высоты, l_rho - плотности.
# Делает:				Выбирает плотность по высоте.
# Возвращает:			Плотность, как это не странно.
def f_rho_dra(altitude, dra):
    
    l_alt = dra['l_alt'].copy()
    l_rho = dra['l_rho'].copy()
    
    if altitude <= 5000:
        ost = math.fmod(altitude, 100)
        if ost>=50:
            altitude_double = altitude + 100 - ost
            index = int(altitude_double/100)
            ret = float(l_rho[index])
            return (ret)
        else:
            altitude_double = altitude - ost
            index = int(altitude_double/100)
            ret = float(l_rho[index])
            return (ret)
    
    if(altitude >=5000 and altitude <= 51000):
        ost = math.fmod(altitude, 500)
        
        if ost < 250:
            altitude_double = altitude - math.fmod(altitude, 500)
            index = 51 + int(altitude_double/500) -3
            return(float(l_rho[index]))
        if ost >= 250:
            altitude_double = altitude - math.fmod(altitude, 500) + 500
            index = 51 + int(altitude_double/500) -3
            return(float(l_rho[index]))
    
    
    if(altitude >= 51000 and altitude < 69500):
        ost = math.fmod(altitude, 1000)
        if(ost < 500):
            altitude_double = altitude - math.fmod(altitude, 1000)
            index = 150 + int(altitude_double/1000) -51
            return(float(l_rho[index]))
        if( ost >= 500):
            altitude_double = altitude - math.fmod(altitude, 1000) + 1000
            index = 150 + int(altitude_double/1000) - 51
            return(float(l_rho[index]))



# Получает:			Флаг взятия решения для программы.
# Делает:			-	Флаг опущен - берёт данные из текстового файла.
# 					-	Флаг поднят - берёт данные из csv-таблицы.
# Возвращает:		Словарь списков: высота <l_alt> и плотность <l_rho>.
def get_table(fcsv=False):
    #Считаем весь файл в список.
    l_alt=[]
    l_rho=[]
    if(fcsv):
        f_out = open('atmosphere_GOST4401_81.csv', 'r')
        for line in f_out:
            line = line.replace('\n', '')
            line = line.replace(':', '')
            lwr=line.split(',')
            l_alt.append(float(lwr[0]))
            l_rho.append(float(lwr[1]))
    else:
        f_out = open('atmosphere_GOST4401_81_data_rho.txt', 'r')
        for line in f_out:
            line = line.replace('\n', '')
            lwr = line.split(',')
            l_alt.append(float(lwr[0]))
            l_rho.append(float(lwr[1]))
    
    d= {
        'l_alt': l_alt,
        'l_rho': l_rho
        }
    return(d)


# Процедура возвращает предельную высоту полёта по таблицам.
def getmax():
    #Считаем весь файл в список.
    l_alt=[]
    l_rho=[]
    f_out = open('atmosphere_GOST4401_81.csv', 'r')
    for line in f_out:
        line = line.replace('\n', '')
        line = line.replace(':', '')
        lwr=line.split(',')
        l_alt.append(float(lwr[0]))
        l_rho.append(float(lwr[1]))
    return(l_alt[-2])



# Экспериментально проверено, что открытие и работа с csv-файлом каждый раз приводит к значительному 
# удлинению времени ответа функции. Работа с открытием текстового файла идёт быстрее в 2.5-3 раза.
# Процедура проверяет наличие текстового файла и разворачивает его при необходимости.
def csv_isfile():
    import os.path
    # Проверка существования.
    flag = os.path.isfile('atmosphere_GOST4401_81_data_rho.txt')
    if(flag == False):
        l_alt=[]
        l_rho=[]
        f_out = open('atmosphere_GOST4401_81.csv', 'r')
        for line in f_out:
            line = line.replace('\n', '')
            line = line.replace(':', '')
            lwr=line.split(',')
            l_alt.append(float(lwr[0]))
            l_rho.append(float(lwr[1]))
        f_out.close()
        
        # Закинуть в считанный файлик данные.
        f_to = open('atmosphere_GOST4401_81_data_rho.txt', 'w')
        counter = -1
        for elem in l_alt:
            counter += 1
            lineto = str(elem) + ',' + str(l_rho[counter]) +'\n'
            f_to.write(lineto)
        f_to.close()
    
    


# Процедура для дополнения плотностей с клавиатуры.
def input_from_keyboard(altitude):    
    
    f_out = open('atmosphere_GOST4401_81.csv', 'a')
    csv.register_dialect('csv_table', delimiter=',', lineterminator='\n')
    fcsv = csv.writer(f_out)
    counter = -1
    while(True):
        counter +=1
        alt = input('alt: ')
        rho = input('rho:')
        lwr = [alt, rho]
        fcsv.writerow(lwr)
    f_out.close()








if __name__ == "__main__":    
    
    flag = True
    l_al = []
    current = 0
    l_al.append(current)
    l_rho = []
    l_rho.append(0)
    while(flag):
        if(current >= 7e4):
            flag = False
        else:
            current += 100
            l_al.append(current)
            l_rho.append(f_rho(current))
    import matplotlib.pyplot as plt
    plt.plot(l_al, l_rho)
    plt.show()
