## Версия 1.
## 
## 	Модуль тригонометрических функций. 
## 
## 
## 


## РАБОТАЮТ ТОЛЬКО МЕТОДЫ sin(teta, flag=False) И cos(teta, flag=False) !!!!


import math


# Получает:				два любых числа
# Делает:				сравнение.
# Возвращает:			модуль разницы между числами.
def delta(x, y):
    if(x==0):
        return(abs(y))
    
    if(y==0):
        return(abs(x))
    
    if(x==y):
        return(0)
    
    if((x>0)and(y>0)):
        #print('x, y ', x, y)
        return(abs(x-y))
    
    if((x<0)and(y>0)):
        return(abs(x) + abs(y))
    
    if(((x<0)and(y<0)) or ((x>0)and(y<0))):
        res = abs(abs(x) - abs(y))
        if(res<0):
            return(res*(-1))
        return(res)

def sin(teta, flag=False):
    if(teta >= 0):
        while(teta >= 360):
            teta -= 360
    if(teta <= 0):
        while(teta<= -360):
            teta += 360
    if(abs(teta) == 0):
        return(0)
    if(flag == False):
        teta_radians = math.radians(teta)
        return math.sin(teta_radians)
    else:
        return math.sin(teta)

def cos(teta, flag=False):
    if(teta >= 0):
        while(teta >= 360):
            teta -= 360
    if(teta <= 0):
        while(teta<= -360):
            teta += 360
    if(abs(teta) == 90):
        return(0)
    if(flag == False):
        teta_radians = math.radians(teta)
        return math.cos(teta_radians)
    else:
        return math.cos(teta)

# Получает:			<y>, type=float, результат применения тригонометрической функции к некоему аргументу. 
# 					<flag>, type=bool, маркер формы возврата решения.
# 					<eps>, type=float, число определяет точность поиска решения.
# Делает:			Вычисляет арксинус заданного числа с заданной точностью.
# Возвращает:		При <flag>==True:
# 						Ответ в долях, пределы: [-1;+1]
# 					При <flag>==False:
# 						Ответ в градусахб пределы: []
def asin(y, flag=False, eps=0.0005):
    if(eps < 0.0005):
        print('ERROR: <trigonometria.py>::asin(%s, flag=%s, eps=%s) - превышена точность выполнения.' %(y, flag, eps))
        return('ERROR: <trigonometria.py>::asin(%s, flag=%s, eps=%s) - превышена точность выполнения.' %(y, flag, eps))
    if(y == 0):
        return(0)
    if(int(y) == 1):
        return(90)
    if(y > 1):
        print("ERROR: y > 1")
        return("ERROR: y > 1")
    elif(y > 0):
        alpha_0 = 0
        alpha_lim = 3.14/2
        while(alpha_0 <= alpha_lim):
            res = delta(sin(alpha_0, True), y)
            if(res<=eps):
                if(flag == True):
                    return(alpha_0)
                return(math.degrees(alpha_0))
            alpha_0 += 0.0005
        print('ERROR: <trigonometria.py>: ошибка в вычислении арксинуса %s' %(y))
        return(None)
    elif(y < 0):
        alpha_0 = 0
        alpha_lim = -3.14/2
        while(alpha_0 >= alpha_lim):
            res = delta(sin(alpha_0, True), y)
            if(res<=eps):
                #print(res)
                if(flag == True):
                    return(alpha_0)
                return(math.degrees(alpha_0))
            alpha_0 -= 0.0005
        print('ERROR: <trigonometria.py>: ошибка в вычислении арксинуса %s' %(y))
        return(None)

def acos(y, flag=False, eps=0.0005):
    if(eps < 0.0005):
        print('ERROR: <trigonometria.py>::acos(%s, flag=%s, eps=%s) - превышена точность выполнения.' %(y, flag, eps))
        return('ERROR: <trigonometria.py>::acos(%s, flag=%s, eps=%s) - превышена точность выполнения.' %(y, flag, eps))
    if(y == 0):
        return(3.14)
    if(abs(y)-eps == 1):
        return(0)
    #if(int(y))
    if(y < -1):
        print('ERROR: <trigonometria.py>::acos(%s, flag=%s, eps=%s) - выход за доступные пределы.' %(y, flag, eps))
        return('ERROR: <trigonometria.py>::acos(%s, flag=%s, eps=%s) - выход за доступные пределы.' %(y, flag, eps))
    if(y > 3.14):
        print("ERROR: y > 1")
        return("ERROR: y > 1")
    elif(y > 0):
        alpha_0 = 0
        alpha_lim = 3.14/2
        while(alpha_0 <= alpha_lim):
            res = delta(cos(alpha_0, True), y)
            if(res<=eps):
                if(flag == True):
                    return(alpha_0)
                return(math.degrees(alpha_0))
            alpha_0 += 0.0005
        print('ERROR: <trigonometria.py>: ошибка в вычислении арккосинуса %s' %(y))
        return(None)
    elif(y < 0):
        alpha_0 = 0
        alpha_lim = -3.14/2
        while(alpha_0 >= alpha_lim):
            res = delta(cos(alpha_0, True), y)
            if(res<=eps):
                #print(res)
                if(flag == True):
                    return(alpha_0)
                return(math.degrees(alpha_0))
            alpha_0 -= 0.0005
        print('ERROR: <trigonometria.py>: ошибка в вычислении арккосинуса %s' %(y))
        return(None)

def tg(teta):
    teta_radians = math.tan(teta)
    return(teta_radians)

def ctg(teta):
    teta_radians = (math.cos(teta) / math.sin(teta))
    return(teta_radians)

def atg(teta):
	si = sin(teta)
	co = cos(teta)
	return(si/co)

#def arctg(teta):
    #teta_radians = math.atan(teta)
    #return(math.degrees(teta_radians))

#def arcctg(teta):
    #teta_radians = math.a


def radians(teta):
    teta = math.radians(teta)
    while(teta>=3.15):
        teta -= 3.14
    return teta


if __name__ == "__main__":
    a = -10
    b = -16
    
    #print(delta(a, b))
    #print(delta(b, a))
    #teta_0 = -360-90
    #print(sin(teta_0))
    #print(cos(teta_0))
    print('asin')
    print(asin(0.3, True))
    print(asin(0.3))
    print(asin(0.3, eps=0.0005))
    print('asin -0.3')
    print(asin(-0.3, True))
    print(asin(-0.3))
    
    import matplotlib.pyplot as plt
    l_sin = []
    l_cos = []
    l_t = []
    t = -360
    while(t <= 360):
        res = cos(t)
        l_cos.append(res)
        res = sin(t)
        l_sin.append(res)
        
        t+= 0.01
        l_t.append(t)
    plt.plot(l_t, l_sin)
    plt.plot(l_t, l_cos)
    plt.legend(['cos', 'sin'])
    plt.grid()
    plt.show()
    
    
    l_asin = []
    l_acos = []
    l_t = []
        
    t = -1
    while(t <= 1):
        res = asin(t)
        l_asin.append(res)
        
        res = acos(t)
        l_acos.append(res)
        
        t+= 0.01
        l_t.append(t)
    #print(l_asin[-1])
    #print(l_t[-1])
    plt.plot(l_t, l_asin)
    plt.plot(l_t, l_acos)
    plt.legend(['asin', 'acos'])
    plt.grid()
    plt.show()
