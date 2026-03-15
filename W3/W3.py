# -*- coding: utf-8 -*-
"""
GF(2^10) Extended Euclidean Algorithm
Irreducible polynomial: m(x) = x^10 + x^3 + 1  (decimal 1033)
Representation: integer where bit i = coefficient of x^i
"""

# ── Hằng số ────────────────────────────────────────────────────────────────
MOD = 0b10000001001  # x^10 + x^3 + 1 = 1033


# ── Tiện ích hiển thị ──────────────────────────────────────────────────────
def poly_str(n: int) -> str:
    """Chuyển số nguyên -> chuỗi đa thức, vd 0b1011 -> 'x^3 + x + 1'."""
    if n == 0:
        return "0"
    terms = []
    i = n.bit_length() - 1
    while i >= 0:
        if n & (1 << i):
            if i == 0:
                terms.append("1")
            elif i == 1:
                terms.append("x")
            else:
                terms.append(f"x^{i}")
        i -= 1
    return " + ".join(terms)


def show(label: str, val: int, indent: int = 4) -> None:
    """In một giá trị ở 3 dạng: decimal / binary / polynomial."""
    pad = " " * indent
    bits = bin(val) if val >= 0 else f"-{bin(-val)}"
    print(f"{pad}{label}")
    print(f"{pad}  dec = {val}")
    print(f"{pad}  bin = {bits}")
    print(f"{pad}  pol = {poly_str(val)}")


# ── Phép toán trên đa thức GF(2) ──────────────────────────────────────────
def degree(n: int) -> int:
    """Bậc của đa thức (bit cao nhất đang bật). -1 nếu n=0."""
    return n.bit_length() - 1


def poly_add(a: int, b: int) -> int:
    """Cộng / trừ đa thức GF(2) = XOR."""
    return a ^ b


def poly_mul(a: int, b: int) -> int:
    """Nhân hai đa thức trên GF(2) (chưa mod)."""
    result = 0
    while b:
        if b & 1:          # nếu bit thấp nhất của b = 1
            result ^= a    # cộng a vào kết quả
        a <<= 1            # dịch a lên 1 bậc (nhân với x)
        b >>= 1            # dịch b xuống
    return result


def divmod_poly(a: int, b: int):
    """
    Chia đa thức a cho b trên GF(2).
    Trả về (q, r) sao cho a = q*b + r, deg(r) < deg(b).
    """
    q = 0
    r = a
    db = degree(b)
    while degree(r) >= db:
        shift = degree(r) - db
        q ^= (1 << shift)   # ghi thêm bit vào thương
        r ^= (b << shift)   # trừ (= XOR) b * x^shift
    return q, r


def poly_mod(a: int, m: int) -> int:
    """Lấy a mod m (chỉ cần phần dư)."""
    _, r = divmod_poly(a, m)
    return r


def gf_mul(a: int, b: int) -> int:
    """Nhân trong GF(2^10): nhân đa thức rồi mod m(x)."""
    return poly_mod(poly_mul(a, b), MOD)


# ── Thuật toán Euclid mở rộng ─────────────────────────────────────────────
def gf_inverse(a: int) -> int:
    """
    Tìm nghịch đảo của a trong GF(2^10) dùng Extended Euclidean.

    Duy trì bất biến:
        r0 = s0 * a_orig  (mod m)
        r1 = s1 * a_orig  (mod m)
    Khi r_i = 1  =>  s_i = a^-1.
    """
    print(f"\n{'='*60}")
    print(f"  Extended Euclidean cho a = {a}")
    print(f"{'='*60}")

    # Khởi tạo: r0 = m(x), r1 = a
    # s0 = 0 (m không liên quan đến a)
    # s1 = 1 (a = 1 * a)
    r0, s0 = MOD, 0
    r1, s1 = a,   1

    print(f"\n[Khởi tạo]")
    show("r0 = m(x)", r0)
    show("r1 = a   ", r1)
    show("s0 = 0   ", s0)
    show("s1 = 1   ", s1)

    step = 0
    while r1 != 0 and r1 != 1:
        step += 1
        print(f"\n{'─'*55}")
        print(f"  Bước {step}")
        print(f"{'─'*55}")

        q, r2 = divmod_poly(r0, r1)

        # s2 = s0 XOR (q*s1 mod m)
        qs1 = poly_mod(poly_mul(q, s1), MOD)
        s2  = poly_add(s0, qs1)

        print()
        show(f"r0", r0)
        show(f"r1", r1)
        show(f"q  = r0 div r1", q)
        show(f"q * r1 (raw, để verify chia)", poly_mul(q, r1))
        show(f"r2 = r0 mod r1", r2)
        show(f"q * s1 (mod m)", qs1)
        show(f"s2 = s0 XOR (q*s1 mod m)", s2)

        # Kiểm tra bất biến: r2 phải bằng s2 * a_orig mod m
        check = gf_mul(s2, a)
        status = "OK" if check == r2 else "FAIL"
        print(f"    [bất biến] s2 * a mod m = {check}  ({status})")

        # Tiến bước
        r0, s0 = r1, s1
        r1, s1 = r2, s2

    if r1 == 0:
        raise ValueError(f"gcd(a, m) != 1 — {a} không có nghịch đảo!")

    # r1 == 1 => s1 là nghịch đảo
    inv = s1
    print(f"\n[Kết thúc] r = 1  =>  a^-1 = s = {inv}")
    show("a^-1", inv)
    return inv


# ── Verify ────────────────────────────────────────────────────────────────
def verify(a: int, inv_a: int) -> None:
    """Kiểm tra a * a^-1 ≡ 1 (mod m)."""
    product = gf_mul(a, inv_a)
    print(f"\n[Verify] {a} * {inv_a} mod m(x)")
    show("product", product)
    status = "PASS ✓" if product == 1 else "FAIL ✗"
    print(f"    => {status}: tích = {product}")


# ── Main ──────────────────────────────────────────────────────────────────
def main():
    print("GF(2^10)  --  m(x) = x^10 + x^3 + 1")
    show("m(x)", MOD, indent=0)

    test_cases = [523, 1015]
    results = {}

    for val in test_cases:
        inv = gf_inverse(val)
        results[val] = inv
        verify(val, inv)

    # ── Tóm tắt kết quả ──────────────────────────────────────────────────
    print("\n" + "="*60)
    print("  KẾT QUẢ CUỐI CÙNG")
    print("="*60)
    for val, inv in results.items():
        print(f"\n  a     = {val:4d}  |  bin = {bin(val)}  |  pol = {poly_str(val)}")
        print(f"  a^-1  = {inv:4d}  |  bin = {bin(inv)}  |  pol = {poly_str(inv)}")
        print(f"  check = {gf_mul(val, inv)}  (phải = 1)")


if __name__ == "__main__":
    main()
