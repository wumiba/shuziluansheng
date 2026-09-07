# -*- coding: utf-8 -*-
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import _fonts

fig, ax = plt.subplots(figsize=(13.6, 5.4))
ax.set_xlim(0, 100); ax.set_ylim(0, 54); ax.axis('off')
fig.canvas.draw()

texts={}
def box(cx,cy,w,h,lines,name,fc='white',ec='#31557a',lw=1.5,fs=8.6,tc='#14181c',bold=False):
    x,y=cx-w/2,cy-h/2
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.02,rounding_size=1.2',fc=fc,ec=ec,lw=lw,zorder=3))
    if isinstance(lines,str): lines=[lines]
    texts[name]=(ax.text(cx,cy,'\n'.join(lines),ha='center',va='center',fontsize=fs,color=tc,zorder=5,linespacing=1.45),cx,cy,w,h)
def arrow(x1,y1,x2,y2,color='#33475b',lw=1.6,ls='-',style='-|>',ms=13):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle=style,mutation_scale=ms,lw=lw,color=color,linestyle=ls,zorder=2,shrinkA=2,shrinkB=2))

# 六阶段主线（单行）
stages = [
 ('S1','三维几何模型','胸腔与心脏三维表面 (STL)','顶点/面片/法向'),
 ('S2','运动机理驱动形变',r'Δ$_i$ = A$_i$·M(t)·$\hat{n}_i$','呼吸/心搏/脉动位移场'),
 ('S3','FMCW 逐点回波','r$_i$(t)、φ$_i$(t)=4πr/λ','77 GHz·快时间采样'),
 ('S4','多散射点接收合成','Mix = Σ g$_i$·e^{jφ_i}','8 通道·帧内积累'),
 ('S5','距离/相位处理','Range-FFT → 相位解调','距离门选取动目标'),
 ('S6','生理量还原与标注','呼吸/心率/PTT/PWV 提取','叠加模型真值标签'),
]
names=[]
w,h = 15.2, 17.0
gap = 1.55
total = 6*w + 5*gap
x0 = (100-total)/2
for i,(nm,t1,t2,t3) in enumerate(stages):
    cx = x0 + w/2 + i*(w+gap)
    fc = ['#e7edf3','#e0f0e7','#ebe3f5','#f3e8ee','#fce9d6','#e3f1e3'][i]
    ec = ['#34506e','#2e6d54','#6b4f9e','#9c5a75','#c2762a','#2f7a3f'][i]
    box(cx, 33, w, h, [t1, t2, t3], nm, fc=fc, ec=ec, fs=8.3, bold=(i==0))
    names.append(nm)
    if i>0:
        axc = x0+w/2+(i-1)*(w+gap)+w/2+gap/2
        arrow(axc-0.2, 33, axc+0.4, 33, color='#445566', lw=2.0, ms=15)
# 循环标注（时间帧）
arrow(x0+w/2, 41.5, x0+w+4.5+0+w+4.5, 41.5, color='#888', ls=':')  # 不实际需要，跳过

# 顶部“逐帧循环”说明
ax.text(50, 52.5, '数字孪生端到端信号仿真链路（每帧：几何刷新 → 回波合成 → 处理；多帧构成时间序列）',
        ha='center', va='center', fontsize=11.5, fontweight='bold', color='#10233a')

# 底部真值同步记录框
box(50, 8.5, 78, 8.5, ['真值标签（由驱动模型与生理参数直接已知，与雷达处理无关）：逐帧位移场 · 呼吸/心率调制 M$_r$/M$_h$ · PTT/PWV · SBP/DBP 设定值',
                       '用于受控评测与机理解释：相位解调/还原结果与模型真值逐帧比对，评估误差并反哺参数标定'], 'TRUTH',
    fc='#fdf1d7', ec='#8a8a3f', fs=8.6)

# 各阶段底部小箭头到真值说明的示意：S2 与 S6 分别拉线

# S2（运动）与 S6 连到真值（虚线）

# 重新用坐标连虚线：S2 cx 与 S6 cx
cx_s2 = x0 + w/2 + 1*(w+gap)
cx_s6 = x0 + w/2 + 5*(w+gap)
arrow(cx_s2, 24.5, cx_s2, 13.0, color='#6b8e23', ls='--', lw=1.3)
arrow(cx_s6, 24.5, cx_s6, 13.0, color='#6b8e23', ls='--', lw=1.3)
ax.text((cx_s2+cx_s6)/2, 10.6, '真值记录与比对', fontsize=7.8, color='#5a6e1f', ha='center')

# 自检溢出
renderer=fig.canvas.get_renderer()
def box_px(cx,cy,w,h):
    x,y=cx-w/2,cy-h/2
    (x0,y0)=ax.transData.transform((x,y)); (x1,y1)=ax.transData.transform((x+w,y+h))
    return min(x0,x1),min(y0,y1),max(x0,x1),max(y0,y1)
print('== overflow check ==')
bad=False
for nm,(t,cx,cy,w,h) in texts.items():
    bb=t.get_window_extent(renderer=renderer)
    bx0,by0,bx1,by1=box_px(cx,cy,w,h)
    wr=bb.width/(bx1-bx0); hr=bb.height/(by1-by0)
    if wr>1 or hr>1:
        print(f'{nm:6s} wr={wr:4.2f} hr={hr:4.2f} OVER!'); bad=True
print('overflow none' if not bad else 'overflow: fix')
print('cx_s2', round(cx_s2,2), 'cx_s6', round(cx_s6,2))

fig.savefig('数字孪生素材_开题报告/图/fig4_回波仿真信号处理链路图.png', dpi=300)
fig.savefig('数字孪生素材_开题报告/图/fig4_回波仿真信号处理链路图.svg')
print('saved fig4')
