#
#	Черноморец Илья
#	Геометрия.
#



# Процедура проводит линию по 2-ум точкам.
# Возвращает коэффициенты в линейном уравнении y = k*x + b
def line2pt(x1, y1, x2, y2):
    k = (y2 - y1)/(x2-x1)
    print('k = ', k)
    
    b = y1 - x1*(y2-y1)/(x2-x1)
    print('b = ', b)


def line2ptin(x1, y1, x2, y2):
    k = (y2 - y1)/(x2-x1)
        
    b = y1 - x1*(y2-y1)/(x2-x1)
    return(k, b)



#	Процедура позволяет построить линейную зависимость по двум точкам и найти на ней нужную точку.
#	Получает:			-	l_x		- массив иксов.
#						-	l_y		- массив игреков.
#						-	x		- координата поиска.
#	Делает:				линеаризует.
#	Возвращает:			значение функции в заданной точке.
def line_get_point(l_x, l_y, x):
	x2 = x
	x1 = x
	counter = -1
	for elem in l_x:
		counter += 1
		if(elem < x):
			x1 = elem
			y1 = l_y[counter]
		if(elem > x):
			x2 = elem
			y2 = l_y[counter]
		if(elem == x):
			return(l_y[counter])
	# Вычисление значения.
	#print(x1, x2, y1)
	k, b = line2ptin(x1, y1, x2, y2)
	return(k*x + b)
    
            
    
if __name__ == "__main__":
    line2pt(0, 0.19, 0.75, 0.19)
    print()
    line2pt(0.75, 0.19, 1.3, 0.39)
    print()
    line2pt(1.3, 0.39, 3.5, 0.28)
    
    print(line_get_point([1, 2], [1, 2], 1.5))
