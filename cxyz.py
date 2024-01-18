#
#	Модуль позволяет создать линеаризованную зависимость коэффициентов и применять её.
#	<cx>, <cy>, <cz>.
#	
#	Известные процедуры:
#	make_cxyz(target)						Создать.
#	get_value(table, xy, flag=True)			Получить значение на псевдокривой.
#	get_table(name)							Получить таблицы из файлов.
            

#	Процедура позволяет создать <cx>, <cy>, <cz> линеаризацией по точкам.
#	Получает:			ничего.
#	Делает:				Создаёт таблицу значений c<x, y, z> по точкам.
#	Возвращает:			1. Создаёт в директории <cxyz> график данного коэффициента.
#						2. Создаёт там же таблицу в формате *.csv.
def make_cxyz(target):
	import libgeometry as libline
	#import matplotlib.pyplot as plt
	import libf2 as libf
	
	# Сохранение данных.
	def save_cx(lx, ly):
		print("==")
		#import matplotlib.pyplot as plt
		import os
		import os.path
		
		# Создание директории и сохранение внутрь файла.
		if(os.path.isdir("cxyz") == False):
			os.mkdir('cxyz')
		
		nameline_alpha = "cxyz_c" + target + '_alpha.txt'
		nameline_c = "cxyz_c" + target + "_c.txt"
		
		libf.L_IN_FILE(nameline_alpha, lx)
		libf.L_IN_FILE(nameline_c, ly)
		libf.files_sort(nameline_alpha, 'cxyz')
		libf.files_sort(nameline_c, 'cxyz')
	# Создание и сохранение картинки графиков.
        #plt.plot(lx, ly)
        #plt.grid()        
        #nameline = "cxyz_c" + target + '.pdf'
        #plt.savefig(nameline, format='pdf', dpi=1000)
        #libf.files_sort(nameline, 'cxyz')
        
        
    
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
		print(l_x)
		print(l_y)
		return(l_x, l_y)
    
    
	print('\nДобро пожаловать в отсек генерации коэффициентов C' + target)
	print('\tВводите свои точки.')

	while(True):
		# 0.
		print('Управление:\tВводите точки, следуете инструкциям.\n')
		
		# 1. Получить список точек.
		lx, ly = inout()

		# 2. Пусть смертные узрят своё безобразное творение!
		#plt.plot(lx, ly)
		#plt.grid()        
		#plt.show()
		
		print('Напоминаем, для выхода СЕЙЧАС нужно нажать <D>.')
		line = input('То ли это? y/n: ')
		if((line == "y")or(line == "Y")):
			save_cx(lx, ly)
			return(lx, ly)
		if((line == 'd') or (line =='D')):
			return(lx, ly)
		else:
			print('Попробуем вбить функцию заново.\n')



#	Процедура позволяет получить значение коэффициента на прямой.
#	Получает:			1. таблицу координат в виде словаря.
#								["x"]	-	коорд. икс.
#								["y"]	-	коорд. игрек.
#						2. координату икс или игрек.
#						3. флаг выбора столбца - 
#							<True> = x,
#							<False> = y.
#	Делает:				Осуществляет выборку данных.
#	Возвращает:			Результат.
def get_value(table, xy, flag=True):
	import libgeometry as libline
	# 1. Считать из таблиц.
	lx = []
	ly = []
	lx = table['x'].copy()
	ly = table['y'].copy()
	#counter =-1
	#for elem in table['x']:
	#    counter +=1
	#    lx.append(elem)
	#    ly.append(table['y'][counter])
	if(type(lx[0]) != float):
		counter = 0
		for elem in lx:
			lx[counter] = float(elem)
			counter += 1
	if(type(ly[0]) != float):
		counter = 0
		for elem in ly:
			ly[counter] = float(elem)
			counter += 1
	if(type(xy) != float):
		xy = float(xy)
	
	# 2. Запустить обсчёт и вернуть значения.
	if(flag==True):
		bufer = libline.line_get_point(lx, ly, xy)
		return(bufer)
	else:
		bufer = libline.line_get_point(ly, lx, xy)
		return(bufer)
    
    

#	Процедура позволяет получить значения таблиц из файлов.
#	Получает:			1. Целевой коэффициент.
#								<cx>	-	получить коэффициент cx
#								<cy>	-	получить cy
#								<cz>	-	получить cz
#	Делает:				Строит таблицу и возвращает её.
#	Возвращает:			таблицу.
def get_table(name):
	import libf2 as libf
	if(name == "cx"):
		lx = libf.L_FROM_FILE('cxyz/cxyz_cx_alpha.txt')
		ly = libf.L_FROM_FILE('cxyz/cxyz_cx_c.txt')
		return(lx, ly)
	if(name == "cy"):
		lx = libf.L_FROM_FILE("cxyz/cxyz_cy_alpha.txt")
		ly = libf.L_FROM_FILE("cxyz/cxyz_cy_c.txt")
		return(lx, ly)
	if(name == 'cz'):
		lx = libf.L_FROM_FILE("cxyz/cxyz_cz_alpha.txt")
		ly = libf.L_FROM_FILE("cxyz/cxyz_cz_c.txt")
		return(lx, ly)




if __name__ == "__main__":
	lx = []
	ly = []
	lx, ly = make_cxyz('x')
	
	#print(get_table("cx"))
