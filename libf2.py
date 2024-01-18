## Версия 2.
## Список методов:

## L_FROM_FILE(filename)
## L_IN_FILE(filename, l_data)

## D_FROM_FILE(filename)
## D_IN_FILE(filename, d_in)
## D_FROM_FILE_FLOAT(filename)
## D_IN_FILE_FLOAT(filename, d_in)
## D_PRINT(d_in)

## files_sort(fname, dirname)
## cleandir(dirname)
## f_remove(fname)

## moving_to_csv		(ПРИМЕР)
## csv_write(filename)
## csv_read(filename)

## input_in_csv()
## 
## 		РИСОВАНИЕ ГРАФИКА И СТРОИТЕЛЬСТВО *.csv!

import os
import shutil
import csv



# Процедура читает данные из файла.
# Получает:				Имя файла для чтения,
# Делает:				Читает файл.
# Возвращает:			Конвертацию строки в число <float>.
def L_FROM_FILE(filename):
    if((filename == None)or(filename == '')):
        return("ERROR: not filename.")
    l_get = get_lines(filename)
    l_ret = convert_lines(l_get)
    return(l_ret)
    
# Процедура читает из файла строки.
# Получает:				имя файла.
# Делает:				Читает и складывает строки в мешок.
# Возвращает:			Список прочитанных строк.
def get_lines(filename):
    f = open(filename, 'r')
    l_str=[]
    for line in f:
        l_str.append(line)
    f.close()
    return(l_str)

# Конвертация.
# Получает:				Список строк.
# Делает:				Конвертирует строку в число с плавающей точкой.
# Возвращает:			Список чисел с плавающей точкой.
def convert_lines(l_str):
    l_convert = []
    for elem in l_str:
        l_elem = elem.split()
        l_convert.append(l_elem[-1])
    return(l_convert)


# Процедура записывает данные в файл столбиком.
# Получает:				1. Имя файла,
# 						2. Список данных.
# Делает:				Записывает данные в файл столбиком.
# возвращает:			<None>
def L_IN_FILE(filename, l_data):
    if((filename == None) or (filename == "")):
        return("ERROR: not filename.")
    f = open(filename, 'w+')
    for elem in l_data:
        f.write(str(elem) + ';')
    f.close()
    return(None)


# Процедура читает словарь данных из файла.
# Получает:				1. Имя файла.
# Делает:				Читает данные, формирует словарь.
# Возвращает:			Словарь данных.
def D_FROM_FILE(filename):
    if((filename == None) or (filename == "")):
        return("ERROR: not filename.")
    f = open(filename, "r")
    l_str = []
    d_out = {}
    for line in f:
        l_str = l_clean(l_str)
        l_str = line.split(':\t')
        l_str[1] = l_str[1].replace('\n', '')
        d_out[str(l_str[0])] = l_str[1]
    f.close()
    return(d_out)

# Процедура записывает словарь в текстовик.
# Получает:				1. Имя файла.
# 						2. словарь данных.
# Делает:				Осуществляет запись данных.
# Возвращает:			Заглушку. 
def D_IN_FILE_FLOAT(filename, d_in):
    if((filename == None) or (filename == "")):
        return("ERROR: not filename.")
    f = open(filename, "w+")
    for key in d_in:
        f.write("%s:\t%s\n" %(key, str(d_in[key])))
    f.close()
    return(None)

def l_clean(l_in):
    del(l_in)
    l_in = []
    return(l_in)

# Процедура читает словарь данных из файла.
# Получает:				1. Имя файла.
# Делает:				Читает данные, формирует словарь.
# Возвращает:			Словарь данных.
def D_FROM_FILE_FLOAT(filename):
	if((filename == None) or (filename == "")):
		return("ERROR: not filename.")
	f = open(filename, "r")
	l_str = []
	d_out = {}
	for line in f:
		l_str = l_clean(l_str)
		l_str = line.split(':\t')
		l_str[1] = float(l_str[1].replace('\n', ''))
		d_out[str(l_str[0])] = l_str[1]
	f.close()
	return(d_out)

# Процедура записывает словарь в текстовик.
# Получает:				1. Имя файла.
# 						2. словарь данных.
# Делает:				Осуществляет запись данных.
# Возвращает:			Заглушку. 
def D_IN_FILE(filename, d_in):
    if((filename == None) or (filename == "")):
        return("ERROR: not filename.")
    f = open(filename, "w+")
    for key in d_in:
        f.write("%s:\t%s\n" %(key, d_in[key]))
    f.close()
    return(None)



# Процедура печатает словарь.
# Получает:				Словарь.
# Делает:				Аккуратно выводит его столбиком.
# Возвращает:			Заглушку.
def D_PRINT(d_in):
    print("D_PRINT")
    if(d_in == None):
        print(None)
        return
    for key in d_in:
        print(key, ":\t", d_in[key])



# Процедура перемещает указанный файл по указанному адресу из текущей директории.
# Получает:					1) имя файла.
# 							2) название целевой директории.
# Делает:					Перемещает.
# Возвращает:				Правду.
def files_sort(fname, dirname):
    file_source = os.getcwd()
    if(os.path.exists(dirname)):
        try:
            os.replace(fname, file_source+'/'+dirname+'/'+fname)
            print('Файл перемещён успешно.')
            return(True)
        except FileNotFoundError:
            print('Возникли проблемы с перемещением файла. Вполне возможно, что это проблема внутри самой библиотеки.')
            return(False)
        #os.replace(fname, file_source+'/'+dirname+'/')
    else:
        os.mkdir(dirname)
        try:
            os.replace(fname, file_source+'/'+dirname+'/'+fname)
            print('Файл перемещён успешно.')
            return(True)
        except FileNotFoundError:
            print('Возникли проблемы с перемещением файла. Вполне возможно, что это проблема внутри самой библиотеки.')
            return(False)
    return(True)


# Процедура очищает заданный каталог путём его удаления и создания нового.
# Получает:				имя каталога.
# Делает:				удаляет и создаёт
# Возвращает:			ничего
def cleandir(dirname):
    #print(os.listdir(dirname))
    ldir = os.listdir(dirname)
    if(ldir == []):
        return()
    
    # Получение полного пути.
    file_source = os.getcwd()
    
    ## Удаление системного файла.
    #hidden_path = file_source+'/'+dirname+'/'+'.directory'
    ##os.remove(hidden_path)
    #f_remove(hidden_path)
    
    # Зачистка каталога.
    for elem in ldir:
        file_source = os.getcwd()
        path = file_source + '/' + dirname + '/' + elem
        os.remove(path)    
    os.rmdir(dirname)
    os.mkdir(dirname)


# Процедура удаляет файл.
# Получает:				имя файла
# Делает:				удаляет
# Возвращает:			-	<True>, если файл был удалён успешно(то есть был найден.)
# 						-	<False>, если файл найден не был.
def f_remove(name):
    file_source = os.getcwd()
    path = file_source + '/' + name
    os.remove(path)    
    return(True)


# Процедура отправляет файл в csv-таблицу. Представлен пример.
# В нём из текстового файла изымаются данные и упаковываются в csv-таблицу.
def moving_to_csv():
    #Считаем весь файл в список.
    import csv
    row0=[]
    row1=[]
    with open("atmosphere_GOST4401_81_data_rho.txt") as f_input:
        for line in f_input:
            line = line.replace('\n', '')
            l_line = line.split(',')
            row0.append(line.split(' ')[0])
            row1.append(line.split(' ')[1])
    
    f_out = open('atmosphere_GOST4401_81.csv', 'w')
    csv.register_dialect('csv_table', delimiter=',', lineterminator='\n')
    fcsv = csv.writer(f_out)
    counter = -1
    while(counter <= len(row0)):
        counter +=1
        lwr = [row0[counter], row1[counter]]
        print(lwr)
        fcsv.writerow(lwr)
    
    f_out.close()

#	Процедура отправляет в *.csv файл заданного имени два столбца данных.
#	Получает:			1. <filename>
#						2. <lx>
#						3. <ly>
#	Делает:				Открывает, записывает.
#	Возвращает:			Ничего.
def csv_write(filename, lx, ly):
    f_out = open(filename, 'w')
    csv.register_dialect('csv_table', delimiter=',', lineterminator='\n')
    fcsv = csv.writer(f_out)
    counter = -1
    while(counter < len(lx)-1):
        counter += 1
        lwr = [lx[counter], ly[counter]]
        fcsv.writerow(lwr)
    f_out.close()

#	Процедура читает *.csv файл заданного имени в два столбца данных.
#	Получает:			1. <filename>
#	Делает:				Открывает, считывает.
#	Возвращает:			1. <lx>
#						2. <ly>
def csv_read(filename):
    import csv
    lx=[]
    ly=[]
    with open(filename) as f_input:
        for line in f_input:
            line = line.replace('\n', '')
            l_line = line.split(',')
            lx.append(line.split(' ')[0])
            ly.append(line.split(' ')[1])
    f_out.close()
    return(lx, ly)

# Процедура позволяет с клавиатуры считать данные и ввести их в csv-таблицу.
# Получает:				ввод с клавиатуры(ручная забивка таблиц)
def input_in_csv():
    f_out = open('atmosphere_GOST4401_81.csv', 'a')
    csv.register_dialect('csv_table', delimiter=',', lineterminator='\n')
    fcsv = csv.writer(f_out)
    counter = -1
    while(True):
        counter +=1
        alt = input('alt: ')
        rho = input('rho:')
        print()
        lwr = [alt, rho]
        fcsv.writerow(lwr)


# Процедура позволяет с клавиатуры считать данные и ввести их в csv-таблицу.
# Получает:				ввод с клавиатуры(ручная забивка таблиц)
# Делает:				вводите, проверяете, забиваете.
# Возвращает:			ничего.
def input_in_csv():
    print('csv-таблица 2-ух столбцовая.')
    tablename = str(input('Введите имя csv-таблицы.'))
    tablename += '.csv'
    f = open(tablename, 'w')
    csv.register_dialect('csv_table', delimiter=',', lineterminator='\n')
    fcsv = csv.writer(f)
    def strin():
        while(True):
            x = input('Coloumn 1 : ')
            y = input('Coloumn 2 : ')
            print()
            strin = x + ',' + y + '\n'
            print(strin)
            print('Чтобы выйти из набора таблицы жмякните x.')
            flag = input('Так будет введено? y/n  >>  ')
            if(flag == 'y'):
                lwr = [x, y]
                return(lwr)
            if(flag == 'x'):
                return(None)
    while(True):
        lwr = strin()
        if(lwr == None):
            f.close()
            return()
        fcsv.writerow(lwr)
        print()
        



#	Процедура позволяет создать <cx>, <cy>, <cz> линеаризацией по точкам.
#	Получает:			ничего.
#	Делает:				Создаёт таблицу значений c<x, y, z> по точкам.
#	Возвращает:			1. Создаёт в директории <cxyz> график данного коэффициента.
#						2. Создаёт там же таблицу в формате *.csv.
def make_cxyz(target):
    import libgeometry as libline
    import matplotlib.pyplot as plt
    import libf2 as libf
    
    # Сохранение данных.
    def save_cx(lx, ly):
        import os
        import os.path
        # Создание директории и сохранение внутрь файла.
        if(os.path.isdir("cxyz") == False):
            os.mkdir('cxyz')
        nameline = "cxyz_c" + target + '.csv'
        libf.csv_write(nameline, lx, ly)
        libf.files_sort(nameline, 'cxyz')
        
        # Создание и сохранение картинки графиков.
        plt.plot(lx, ly)
        plt.grid()        
        nameline = "cxyz_c" + target + '.pdf'
        plt.savefig(nameline, format='pdf', dpi=1000)
        libf.files_sort(nameline, 'cxyz')
        
        
    
    # Ввод данных в диалаговом режиме и их редактирование.
    def inout():
        l_x = []
        l_y = []
        while(True):
            line = input('Дальше? Y/N ')
            if((line == "N")or(line == "n")):
                break
            if((line == 'y')or(line == 'Y')):
                line = input('Введите x точки: ')
                if(line.isdigit() != True):
                    while(line.isdigit() == False):
                        print('Вводите коряво. ')
                        line = input('Введите x точки: ')
                
                x = float(line)
                
                line = input('Введите y точки: ')
                if(line.isdigit() != True):
                    while(line.isdigit() == False):
                        print('Вводите коряво. ')
                        line = input('Введите y точки: ')
                
                y = float(line)
                
                line_to = 'Правильно ли: '+ str(x) + ' ' + str(y) + ' ? Y/N: '
                flag = bool(input(line_to))
                if(flag):
                    l_x.append(x)
                    l_y.append(y)
                else:
                    print('Тогда на следующей точке попробуйте вбить эту же.')
        return(l_x, l_y)
    
    
    print('\nДобро пожаловать в отсек генерации коэффициентов C' + target)
    print('\tВводите свои точки.')
    
    while(True):
        # 0.
        print('Управление:\tВводите точки, следуете инструкциям.\n')
        # 1. Получить список точек.
        lx, ly = inout()
        
        # 2. Пусть смертные узрят своё безобразное творение.
        plt.plot(lx, ly)
        plt.grid()        
        plt.show()
        print('Напоминаем, для выхода СЕЙЧАС нужно нажать <D>.')
        line = input('То ли это? y/n: ')
        if((line == "y")or(line == "Y")):
            save_cx(lx, ly)
            return(lx, ly)
        if((line == 'd') or (line =='D')):
            return(lx, ly)
        else:
            print('Попробуем вбить функцию заново.\n')




if __name__ == "__main__":
    input_in_csv()
    
    
    
    #filename = "ВВЕДИТЕ ВАШЕ ИМЯ.txt"
    #print(L_FROM_FILE(filename))
    #l_data = L_FROM_FILE(filename)
    
    #filename = "ВВЕДИТЕ ВАШЕ ИМЯ.txt"
    #L_IN_FILE(filename, l_data)
    #f = open(filename, 'r')
    #for line in f:
        #print(line)
    
    #filename = "ВВЕДИТЕ ВАШЕ ИМЯ.txt"
    #d_out = D_FROM_FILE(filename)
    #print(d_out)
    
    #filename = "ВВЕДИТЕ ВАШЕ ИМЯ.txt"
    #D_IN_FILE(filename, d_out)
    
    #D_PRINT(d_out)
    
    
