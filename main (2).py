import random
import pygame
import pygameMenu
import sys
import time

from pygameMenu.locals import * 

import random
FINISH = False

def text1(word,x,y):
    font = pygame.font.SysFont(None, 50)
    text = font.render("{}".format(word), True, (255,255,255))
    a.blit(text,(x,y))
def text2(word,x,y):
    font = pygame.font.SysFont(None, 50)
    text = font.render("{}".format(word), True, (255, 0, 0))
    a.blit(text,(x,y))
def text3(word,x,y):
    font = pygame.font.SysFont(None, 50)
    text = font.render("{}".format(word), True, (255,255,255))
    a.blit(text,(x,y))
def text4(word,x,y):
    font = pygame.font.SysFont(None, 50)
    text = font.render("{}".format(word), True, (0, 0, 0))
    a.blit(text,(x,y))    
def inpt():
    word=""
    pygame.display.flip()
    done = True
    while done:
        a.fill((0,0,0))
        text1("Please enter your name : ",300,400)
        for event in pygame.event.get():
            if event.type==pygame.QUIT:
                pygame.quit()
                quit()
            if event.type == pygame.KEYDOWN:
                if event.key in range(97, 123):
                    word+=str(chr(event.key))
                if event.key == pygame.K_BACKSPACE:
                    word=word[:int(len(word)-1)]
                if event.key == pygame.K_RETURN:
                    done=False
                if len(word)==1:
                    word=word.capitalize()
                #events...
            text1(word,720,400)
            pygame.display.update()
    return word
    
def solo():
   
    
    width=1600
    length=900
    
    
    a=pygame.display.set_mode((0,0),pygame.FULLSCREEN)
    pygame.display.set_caption('Just Game')
    
    score=0
    word=inpt()
    
    x,y,r=width//2,length//2,50
    vx=random.choice([-5,5])
    vy=random.choice([-2,2])
    ltx,lty=0,0
    rwidth,rheight=20,150
    
    helpbon1=0
    hv1=0    
    
    flbon=-1
    
    countercol=0
    counter_col_circ=0
    
    start_timer=time.time()
    
    clock = pygame.time.Clock()
    ska=0
    br = False
    while not(br):
        ska+=1
        if countercol<10:
            a.fill((0,0,0))
            color=(0,0,0)
            chcolor=(255,255,255)
            text3(word,0,0)
            text3(score,0,35)           
        if 20>countercol>=10:
            a.fill((255,255,255))
            color=(255,255,255)
            chcolor=(0,0,0)
            text4(word,0,0)
            text4(score,0,35)           
        if countercol==20:
            countercol=0
            color=(0,0,0)
            chcolor=(255,255,255)        
       
        
        rectangle=pygame.draw.rect(a, (chcolor), (ltx,lty,rwidth,rheight))
    
        curtime=time.time()
        ticks=(curtime-start_timer)
    
        if flbon==1:
            bonrect1=pygame.draw.rect(a, (199, 125, 199), (0,yr,20,20))
            if rectangle.colliderect(bonrect1):
                helpbon1+=1
                hv1=1
            if hv1==1:
                ancdrect=pygame.draw.rect(a, (color), (0,yr,20,20))
                rectangle=pygame.draw.rect(a, (chcolor), (ltx,lty,rwidth,rheight))
        
            if helpbon1>0:
                score+=10
                helpbon1=-123972481759            
        millis = int(str(ticks-int(ticks))[2:4])
    
        secondsw=int(ticks)%60
        secondsb=int(ticks)
    
        minutes=secondsb//60
        out='{minutes:02d}:{seconds:02d}:{millis}'.format(minutes=minutes,seconds=secondsw,millis=millis)
        text2(out, 750,0)                
        if ska==1:
            circcol=(128,128,128)
        if counter_col_circ==7:
            circcol=random.choice([(128, 0, 255),(0, 255, 0),(0, 128, 255)])
            counter_col_circ=0
        circle=pygame.draw.circle(a,circcol, (x, y), r)
            
        
        for event in pygame.event.get():
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_ESCAPE:            
                    br = True
                    ticks=0
        keys = pygame.key.get_pressed()
        if keys[pygame.K_DOWN]:
            lty+=10
            if lty>length-rheight:
                lty=length-rheight            
        if keys[pygame.K_UP]:
            lty-=10
            if lty<0:
                lty = 0            
        x+=vx
        y+=vy
        if y+r>=length:
            vy=-vy
        if y-r<=0:
            vy=-vy
        if x+r>=width:
            flbon=random.randint(1,3)
            vx=-vx
            yr=random.randint(70,875)
            help_rect=pygame.draw.rect(a, (color), (0,yr,20,20))
            while rectangle.colliderect(help_rect)!=False:
                yr=random.randint(70,875)
                help_rect=pygame.draw.rect(a, (color), (0,yr,20,20))
                
            helpbon1=0
            hv1=0
        if x<=r+rwidth and (y-r<=lty+rheight and lty<=y+r):
            flbon=0
            countercol+=1
            counter_col_circ+=1
            if abs(vx)<20:
                vx=abs(vx)+1
            else:
                vx=-vx   
            if lty+rheight/6>=y+r>=lty or lty+rheight*5/6<=y-r<=lty+rheight:
                score+=10
            elif lty+rheight/6<y+r<=lty+rheight*2/6 or lty+rheight*5/6>y-r>=lty+rheight*4/6:
                score+=7.5
            else:
                score+=5
        if x-r<=0:
            a.fill((0,0,0))
            text1("Game Over! Don't worry, you will be lucky next time!",325,400)
            text1("Your score : {}".format(score),650 ,500)
            text1('Total time: {}'.format(out),600,700)
            text1("Press Esc to go to Main Menu",525,800)
            pygame.display.update()  
            while not(br):
                for event in pygame.event.get():
                    if event.type == pygame.KEYDOWN:
                            if event.key == pygame.K_ESCAPE:                                          
                                br = True
                                FINISH = True
        clock.tick(100)
        pygame.display.update()

def doword():
    word1=""
    word2=""
    pygame.display.flip()
    done = True
    w_n = 1
      
    while done:
        for event in pygame.event.get():
            if event.type==pygame.QUIT:
                pygame.quit()
                quit()
            if event.type == pygame.KEYDOWN:
                if event.key in range(97, 123):
                    if w_n == 1:
                        word1+=str(chr(event.key))  
                    else:
                        word2+=str(chr(event.key))  
                if len(word1)==1:
                    word1=word1.capitalize()                
                if len(word2)==1:
                    word2=word2.capitalize()                
                if event.key == pygame.K_BACKSPACE:
                    if w_n == 1:
                        word1=word1[:int(len(word1)-1)]
                    else:
                        word2=word2[:int(len(word2)-1)]
                if event.key == pygame.K_RETURN:
                    w_n+=1
                    if w_n>2:
                        done=False
                #events...
            a.fill((0,0,0))
            text1("Right player, please, enter your name : ",300,650)
            text1("Left player, please, enter your name : ",300,400)    
            text1(word1,920,400)
            text1(word2,945,650)
            pygame.display.update()
    return word1, word2         

flag=0
flag1=0

def inpt1():
    word1=""
    word2=""
    a.fill((0,0,0))
    pygame.display.flip() 
    word1, word2 = doword()
    pygame.display.update()
    return word1,word2

def duo():
    
    width=1600
    length=900
    
    a=pygame.display.set_mode((0,0), pygame.FULLSCREEN)
    pygame.display.set_caption('Just Game')
    score1=0
    score2=0
    words=inpt1()
    word1=words[0]
    word2=words[1]
    x,y,r=width//2,length//2,50
    vx=random.choice([-5,5])
    vy=random.choice([-2,2])
    
    
   
    rwidth1,rheight1=20,150
    rwidth2,rheight2=20,150
    ltx1,lty1=0,0
    ltx2,lty2=width-rwidth2,0  
    
    flbon1=-1
    flbon2=-1
    helpbon1=0
    helpbon2=0    
    hv1=0
    hv2=0
    
    countercol=0
    
    start_timer=time.time()
    counter_col_circ=0
    
    clock = pygame.time.Clock()
    br = False
    ska=0    
    while not(br):
        ska+=1
        if countercol<10:
            a.fill((0,0,0))
            color=(0,0,0)
            chcolor=(255,255,255)
            text3(word1,0,0)
            text3(word2,width-(17*len(word2)),0)
            text3(score1,0,35)
            text3(score2,width-(17*len(str(score2))),35)               
        if 20>countercol>=10:
            a.fill((255,255,255))
            color=(255,255,255)
            chcolor=(0,0,0)
            text4(word1,0,0)
            text4(word2,width-(17*len(word2)),0)
            text4(score1,0,35)
            text4(score2,width-(17*len(str(score2))),35)               
        if countercol==20:
            countercol=0
            color=(0,0,0)
            chcolor=(255,255,255)
                      
        
        curtime=time.time()
        ticks=(curtime-start_timer)
    

        millis = int(str(ticks-int(ticks))[2:4])
        
        secondsw=int(ticks)%60
        secondsb=int(ticks)
        
        minutes=secondsb//60
        out='{minutes:02d}:{seconds:02d}:{millis}'.format(minutes=minutes,seconds=secondsw,millis=millis)
        text2(out, 750,0)        
        
        
        
        if countercol<10:
            rectangle1=pygame.draw.rect(a, (255, 255, 255), (ltx1,lty1,rwidth1,rheight1))
            rectangle2=pygame.draw.rect(a, (255, 255, 255), (ltx2,lty2,rwidth2,rheight2))           
        if 20>countercol>=10:
            rectangle1=pygame.draw.rect(a, (0, 0, 0), (ltx1,lty1,rwidth1,rheight1))
            rectangle2=pygame.draw.rect(a, (0, 0, 0), (ltx2,lty2,rwidth2,rheight2))            
        if countercol==20:
            countercol=0        
       
        if ska==1:
            circcol=(128,128,128)
        if counter_col_circ==7:
            circcol=random.choice([(128, 0, 255),(0, 255, 0),(0, 128, 255)])
            counter_col_circ=0
        circle=pygame.draw.circle(a,circcol, (x, y), r)
        if flbon1==1:
            if countercol<10:
                bonrect1=pygame.draw.rect(a, (251, 255, 0), (0,y1r,20,20))
            else:
                bonrect1=pygame.draw.rect(a, (34, 168, 70), (0,y1r,20,20))            
            if rectangle1.colliderect(bonrect1):
                helpbon1+=1
                hv1=1
            if hv1==1:
                ancdrect1=pygame.draw.rect(a, (color), (0,y1r,20,20))
                rectangle1=pygame.draw.rect(a, (chcolor), (ltx1,lty1,rwidth1,rheight1))
              
            if helpbon1>0:
                score1+=10
                helpbon1=-123972481759                    
        
        if flbon2==1:
            if countercol<10:
                bonrect2=pygame.draw.rect(a, (251, 255, 0), (1580,y2r,20,20))
            else:
                bonrect2=pygame.draw.rect(a, (34, 168, 70), (1580,y2r,20,20))
            if rectangle2.colliderect(bonrect2):
                helpbon2+=1
                hv2=1
            if hv2==1:
                ancdrect2=pygame.draw.rect(a, (color), (1575,y2r,20,20))
                rectangle2=pygame.draw.rect(a, (chcolor), (ltx2,lty2,rwidth2,rheight2))
            
            if helpbon2>0:
                score2+=10
                helpbon2=-123972481759
        
        pygame.display.update() 
        for event in pygame.event.get():
            if event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_ESCAPE:            
                        br = True
        keys = pygame.key.get_pressed()
        if keys[pygame.K_w]:
            lty1-=10
            if lty1<0:
                lty1 = 0
        if keys[pygame.K_s]:
            lty1+=10
            if lty1>length-rheight1:
                lty1=length-rheight1
        if keys[pygame.K_UP]:
            lty2-=10
            if lty2<=0:
                lty2=0
        if keys[pygame.K_DOWN]:
            lty2+=10
            if lty2>=length-rheight2:
                lty2=length-rheight2
        x+=vx
        y+=vy
        if y+r>=length:
            vy=-vy
        if y-r<=0:
            vy=-vy
        if x>=width-r-rwidth2 and (lty2<y-r<lty2+rheight2 or lty2<y+r<lty2+rheight2):
            counter_col_circ+=1
            countercol+=1
            flbon1=random.randint(1,3)
            flbon2=-1
            helpbon1=0
            hv1=0
            y1r=random.randint(70,875)
            while lty1+rheight1>=y1r>=lty1-20:
                y1r=random.randint(70,875)            
            if abs(vx)<20:        
                vx=-abs(vx)-1
            else:
                vx=-abs(vx)
            if lty2+rheight2/6>=y+r>=lty2 or lty2+rheight2*5/6<=y-r<=lty2+rheight2:
                score2+=10
            elif lty2+rheight2/6<y+r<=lty2+rheight2*2/6 or lty+rheight2*5/6>y-r>=lty2+rheight2*4/6:
                score2+=7.5
            else:
                score2+=5  
        elif x<=r+rwidth1 and (lty1<y-r<lty1+rheight1 or lty1<y+r<lty1+rheight1):
            flbon2=random.randint(1,3)
            flbon1=-1
            helpbon2=0
            hv2=0
            y2r=random.randint(70,875)
            while lty2+rheight2>=y2r>=lty2-20:
                y2r=random.randint(70,875)
            if abs(vx)<20:    
                vx=abs(vx)+1
            else:
                vx=abs(vx)
            if lty1+rheight1/6>=y+r>=lty1 or lty1+rheight1*5/6<=y-r<=lty1+rheight1:
                score1+=10
            elif lty1+rheight1/6<y+r<=lty1+rheight1*2/6 or lty1+rheight1*5/6>y-r>=lty1+rheight1*4/6:
                score1+=7.5
            else:
                score1+=5    
        if x-r<=0:
            a.fill((0,0,0))
            score1=score1*0.9
            if score1>score2:
                text1('{} is a winner! Congrats!'.format(word1),550,400)
                text1('{} score : {}'.format(word1,score1),350,600)
                text1('{} score : {}'.format(word2,score2),1000,600)
                text1('Total time: {}'.format(out),600,700)
                text1("Press Esc to go to Main Menu",525,800)            
            if score2>score1:
                text1('{} is a winner! Congrats!'.format(word2),550,400)
                text1('{} score : {}'.format(word1,score1),350,600)
                text1('{} score : {}'.format(word2,score2),1000,600)
                text1('Total time: {}'.format(out),600,700)
                text1("Press Esc to go to Main Menu",525,800)                 
            if score1==score2:
                text1("It's a draw!",550,400)
                text1('{} score : {}'.format(word1,score1),350,600)
                text1('{} score : {}'.format(word2,score2),1000,600)
                text1('Total time: {}'.format(out),600,700)
                text1("Press Esc to go to Main Menu",525,800)                       
            pygame.display.update()  
            while not(br):
                for event in pygame.event.get():
                    if event.type == pygame.KEYDOWN:
                            if event.key == pygame.K_ESCAPE:                    
                                br = True
                                FINISH=True
        if x+r>=width:
            a.fill((0,0,0))
            score2=score2*0.9
            if score1>score2:
                text1('{} is a winner! Congrats!'.format(word1),550,400)
                text1('{} score : {}'.format(word1,score1),350,600)
                text1('{} score : {}'.format(word2,score2),1000,600)
                text1('Total time: {}'.format(out),600,700)
                text1("Press Esc to go to Main Menu",525,800)            
            if score2>score1:
                text1('{} is a winner! Congrats!'.format(word2),550,400)
                text1('{} score : {}'.format(word1,score1),350,600)
                text1('{} score : {}'.format(word2,score2),1000,600)
                text1('Total time: {}'.format(out),600,700)
                text1("Press Esc to go to Main Menu",525,800)                 
            if score1==score2:
                text1("It's a draw!",550,400)
                text1('{} score : {}'.format(word1,score1),350,600)
                text1('{} score : {}'.format(word2,score2),1000,600)
                text1('Total time: {}'.format(out),600,700)
                text1("Press Enter to go to Main Menu",525,800)              
            pygame.display.update()  
            while not(br):
                for event in pygame.event.get():
                    if event.type == pygame.KEYDOWN:
                            if event.key == pygame.K_ESCAPE:                    
                                br = True
                                FINISH=True                
        clock.tick(100)  

width=1600
length= 900
a=pygame.display.set_mode((0,0),pygame.FULLSCREEN)
pygame.display.set_caption('Just Game')

x,y,r=300,400,50
vx=5
vy=2
ltx,lty=0,0
rwidth,rheight=20,150
pygame.init()
clock = pygame.time.Clock()

def mainmenu_background():
    """
    Background color of the main menu, on this function user can plot
    images, play sounds, etc.
    """
    a.fill((40, 40, 40))

HELP = ['Press Enter to enable/disable Menu',
        'If you are playing solo : press UP/DOWN to move slider',
        'If you are playing duo: left player use W/S, right - UP/DOWN to move slider',
        'If you want to exit: choose "Exit" and press Enter in Main Menu']

help_menu = pygameMenu.TextMenu(a,
                                dopause=False,
                                font=pygameMenu.font.FONT_FRANCHISE,
                                menu_color=(0, 0, 0),  # Background color
                                menu_color_title=(120, 45, 30),
                                #onclose=PYGAME_MENU_CLOSE,  # Pressing ESC button does nothing
                                title='Help',
                                window_height=900,
                                window_width=1600
                                )
help_menu.add_option('Return to Menu', pygameMenu.events.BACK)
for m in HELP:
    help_menu.add_line(m)

menu = pygameMenu.Menu(a,
                           bgfun=mainmenu_background,
                           enabled=True,
                           font=pygameMenu.font.FONT_NEVIS,
                           menu_alpha=100,
                           #onclose=PYGAME_MENU_CLOSE,
                           title='Main Menu',
                           title_offsety=8,
                           window_height=length,
                           window_width=width)
menu.add_option("Solo", solo)                                             
menu.add_option("Duo", duo) 
menu.add_option("Help", help_menu) 
menu.add_option('Exit',  pygameMenu.events.EXIT) 



clock.tick(160)



while True:
    events = pygame.event.get()
    
    menu.mainloop(events)
    
    
pygame.display.quit()
pygame.quit() 