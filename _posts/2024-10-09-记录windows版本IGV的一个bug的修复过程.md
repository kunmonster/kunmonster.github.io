---
layout: post
title: 记录windows版本IGV的一个bug的修复过程
date: 2024-10-09 15:11
tags: [bioinformatics,igv]
toc: true
math: true
---



## 问题
---

该bug在[issue #1512](https://github.com/igvteam/igv/issues/1512)中第一次被提及，具体体现为在windows平台下打开IGV，显示IGV的命令行窗口但是不显示IGV主窗口，并且在状态栏有IGV的图标显示。

本人在使用时遇到上述同样的情况，并且经过一段时间定位，发现问题出现在窗口坐标上。问题具体成因如下：

- IGV通过保存一个偏好文件以达到窗口在下次打开的时候重新显示在上次关闭的同一位置
- 当有多块屏幕时，如果上次关闭时IGV所在位置的屏幕没有被连接（特别是作为主屏幕左边的副屏时）那么此时将会出现上述情况
- 根本原因就是指定窗口大小和位置的(x,y,width,height)在上次退出时保存为
  $$ x+width <0  $$
  或者
  $$ y+height < 0 $$
导致在当前屏幕上不能显示出该窗口。

## 解决
---

由于作者使用Mac OS，无法复现该问题，所以我计划修复此bug。


最初我计划给 $$ x+width $$ 和 $$ y+height $$设定一个最小值，如果小于该最小值那么直接将坐标强行拉回(0,0,1150,800)，但是此种做法将会改变原本的保留用户位置偏好的功能。在与作者讨论并且其建议采用以下修复方法

- 获取所有显示设备窗口坐标并存于数组中
- 读取用户偏好坐标，并且判断左上角(x,y)是否在某一个显示设备内，如果不在，那么在主屏幕显示，并且使用默认窗口大小，如果在，那么判断用户偏好中整个窗口是否在该屏幕内，如果不在那么调整窗口大小，显示在该屏幕中

具体实现见 commit 1d92ce47fc6e319222ae942f651fa14a15833e35 in [issue #1590](https://github.com/igvteam/igv/pull/1590)

```java 
GraphicsEnvironment graphEnv = GraphicsEnvironment.getLocalGraphicsEnvironment();  
GraphicsDevice[] graphDev = graphEnv.getScreenDevices();  
Rectangle[] boundsArr = new Rectangle[graphDev.length];  
  
for (int i=0;i<graphDev.length;++i) {  
    GraphicsConfiguration curCon = graphDev[i].getDefaultConfiguration();  
    boundsArr[i] = curCon.getBounds();  
}  
  
//set a flag which indicates if the user preference is empty, or if the (x,y) in the user preference is not contained in any screen  
//default is empty or not contained  
boolean isNullOrNotContained = true;  
  
if(applicationBounds != null){  
    //Iterate over each screen value to find if there is currently a screen that can contain these values.  
    int userX = applicationBounds.x;  
    int userY = applicationBounds.y;  
    double userMaxX = applicationBounds.getMaxX();  
    double userMaxY = applicationBounds.getMaxY();  
    for(Rectangle curScreen : boundsArr){  
        if(curScreen.contains(userX,userY)){  
            isNullOrNotContained = false;  
            if( userMaxX >= curScreen.getMaxX() || userMaxY >= curScreen.getMaxY()){  
                applicationBounds = new Rectangle(curScreen.x,curScreen.y,Math.min(1150,curScreen.width),Math.min(800,curScreen.height));  
            }  
            break;  
        }  
    }  
}  
if(isNullOrNotContained){  
    // user's preference is null or the (x,y) in user's preference is not contained in any screen  
    // set the application to the main screen    
    applicationBounds = new Rectangle(0, 0, Math.min(1150,screenBounds.width),
    Math.min(800,screenBounds.height));  
}
```


该PR已经被作者从fix_issue1512分支合并至main分支
