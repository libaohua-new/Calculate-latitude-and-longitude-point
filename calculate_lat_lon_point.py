import numpy as np
from pyproj import Geod, Transformer

# --------------------------------------------------
# 1. 度分秒 → 十进制度
# --------------------------------------------------
def dms_to_decimal(d: int, m: int, s: float) -> float:
    return d + m / 60 + s / 3600

# 已知 A、B 的经纬度（WGS84）
A_lon = dms_to_decimal(105, 29, 52.3439)   # 东经
A_lat = dms_to_decimal(39, 17, 53.4608)    # 北纬
B_lon = dms_to_decimal(105, 29, 52.3031)   # 东经
B_lat = dms_to_decimal(39, 17, 53.9129)    # 北纬

# C 到 A、B 的距离（米）
dA = 16.91
dB = 9.0

# --------------------------------------------------
# 2. 建立高精度测地对象与局部投影
# --------------------------------------------------
geod = Geod(ellps='WGS84')                 # 椭球面计算
local_tm = Transformer.from_crs(
    {"proj": 'latlong', "datum": 'WGS84'},
    {"proj": 'tmerc', "lat_0": A_lat, "lon_0": A_lon, "datum": 'WGS84'},
    always_xy=True
)

# --------------------------------------------------
# 3. 经纬度 ↔ 局部平面坐标（精确）
# --------------------------------------------------
def latlon_to_local(lon: float, lat: float) -> tuple[float, float]:
    """返回以 A 为原点的局部 TM 坐标（东、北，单位 m）"""
    e, n = local_tm.transform(lon, lat)
    e0, n0 = local_tm.transform(A_lon, A_lat)
    return e - e0, n - n0

def local_to_latlon(e: float, n: float) -> tuple[float, float]:
    """局部 TM 坐标 → 经纬度（度）"""
    e0, n0 = local_tm.transform(A_lon, A_lat)
    lon, lat = local_tm.transform(e + e0, n + n0, direction='INVERSE')
    return lat, lon

# A、B 的局部平面坐标
xA, yA = 0.0, 0.0
xB, yB = latlon_to_local(B_lon, B_lat)

# --------------------------------------------------
# 4. 两圆交点（平面直角坐标系）
# --------------------------------------------------
def circle_intersection(x1, y1, r1, x2, y2, r2):
    d = np.hypot(x2 - x1, y2 - y1)
    if d > r1 + r2 or d < abs(r1 - r2):
        raise ValueError("无解：两圆不相交")
    a = (r1 ** 2 - r2 ** 2 + d ** 2) / (2 * d)
    h = np.sqrt(r1 ** 2 - a ** 2)
    xm = x1 + a * (x2 - x1) / d
    ym = y1 + a * (y2 - y1) / d
    dx = h * (y2 - y1) / d
    dy = h * (x1 - x2) / d
    return (xm + dx, ym + dy), (xm - dx, ym - dy)

try:
    C1_local, C2_local = circle_intersection(xA, yA, dA, xB, yB, dB)
except ValueError as e:
    print(e)
    exit()

# --------------------------------------------------
# 5. 转回经纬度
# --------------------------------------------------
C1_lat, C1_lon = local_to_latlon(*C1_local)
C2_lat, C2_lon = local_to_latlon(*C2_local)

# --------------------------------------------------
# 6. 结果输出
# --------------------------------------------------
def decimal_to_dms(deg):
    d = int(deg)
    m = int((deg - d) * 60)
    s = ((deg - d) * 60 - m) * 60
    return d, m, s

print("两个可能的解（十进制度）：")
print(f"解1：({C1_lat:.8f}, {C1_lon:.8f})")
print(f"解2：({C2_lat:.8f}, {C2_lon:.8f})")

print("\n两个可能的解（度分秒）：")
print("解1：纬度 {:02d}°{:02d}'{:08.5f}\"，经度 {:02d}°{:02d}'{:08.5f}\"".format(
    *decimal_to_dms(C1_lat), *decimal_to_dms(C1_lon)))
print("解2：纬度 {:02d}°{:02d}'{:08.5f}\"，经度 {:02d}°{:02d}'{:08.5f}\"".format(
    *decimal_to_dms(C2_lat), *decimal_to_dms(C2_lon)))

# --------------------------------------------------
# 7. （可选）验证距离
# --------------------------------------------------
for idx, (lat, lon) in enumerate([(C1_lat, C1_lon), (C2_lat, C2_lon)], 1):
    _, _, dist_A = geod.inv(A_lon, A_lat, lon, lat)
    _, _, dist_B = geod.inv(B_lon, B_lat, lon, lat)
    print(f"\n验证 —— 解{idx} 到 A 的距离：{dist_A:.3f} m，到 B 的距离：{dist_B:.3f} m")
