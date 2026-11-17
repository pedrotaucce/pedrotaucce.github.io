#!/usr/bin/env python3
"""Calcula o enquadramento dos cartões de People a partir da posição do rosto.

Para cada foto em assets/img/people/ o script procura o maior rosto (detector YuNet, do OpenCV)
e calcula o `object-position` CSS para o recorte quadrado do cartão ficar centrado no rosto.
O resultado vai para _data/people_photo_positions.yml, que a página People usa quando a pessoa
não tem `photo_position` escrito à mão em _data/people.yml.

O script NUNCA deve derrubar o deploy: se algo faltar (OpenCV, modelo) ou der erro, ele avisa,
não altera nada e termina com sucesso. Fotos sem rosto detectado simplesmente ficam de fora
(a página usa a posição padrão).

Uso (precisa de opencv-python-headless): python3 scripts/photo_positions.py
"""
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MODEL = os.path.join(HERE, "models", "face_detection_yunet_2023mar.onnx")
PHOTOS = os.path.join(ROOT, "assets", "img", "people")
OUT = os.path.join(ROOT, "_data", "people_photo_positions.yml")
EXTENSIONS = (".jpg", ".jpeg", ".png", ".webp")
MIN_SCORE = 0.6          # confiança mínima do detector
FACE_HEIGHT_FRACTION = 0.40  # em fotos em pé, o rosto fica a 40% da altura do quadrado


def position_for(width, height, face):
    """object-position (x%, y%) para um recorte quadrado (object-fit: cover) centrado no rosto."""
    x, y, w, h = face[:4]
    cx, cy = x + w / 2, y + h / 2
    side = min(width, height)
    if width == height:
        return None  # foto quadrada: o enquadramento não importa
    px = 50.0
    py = 50.0
    if width > side:    # foto deitada: só o deslocamento horizontal conta
        px = (cx - side / 2) / (width - side) * 100
    if height > side:   # foto em pé: só o deslocamento vertical conta
        py = (cy - FACE_HEIGHT_FRACTION * side) / (height - side) * 100
    clamp = lambda v: int(round(min(100.0, max(0.0, v))))
    return clamp(px), clamp(py)


def load_previous():
    """Lê os valores já gravados (para mantê-los se uma foto der erro)."""
    previous = {}
    try:
        with open(OUT, encoding="utf8") as handle:
            for line in handle:
                line = line.strip()
                if line.startswith('"') and '": "' in line:
                    name, value = line.split('": "', 1)
                    previous[name[1:]] = value.rstrip('"')
    except OSError:
        pass
    return previous


def main():
    try:
        import cv2
    except ImportError:
        print("photo_positions: OpenCV não está instalado; nada foi alterado.")
        return 0
    if not os.path.exists(MODEL) or not hasattr(cv2, "FaceDetectorYN"):
        print("photo_positions: modelo ou detector indisponível; nada foi alterado.")
        return 0

    previous = load_previous()
    results = {}
    for path in sorted(glob.glob(os.path.join(PHOTOS, "*"))):
        name = os.path.basename(path)
        if not name.lower().endswith(EXTENSIONS):
            continue
        try:
            image = cv2.imread(path)  # respeita a orientação EXIF, como o navegador
            if image is None:
                print(f"photo_positions: não consegui abrir {name}; mantido o valor anterior, se houver.")
                if name in previous:
                    results[name] = previous[name]
                continue
            height, width = image.shape[:2]
            detector = cv2.FaceDetectorYN.create(MODEL, "", (width, height), MIN_SCORE, 0.3, 5000)
            _, faces = detector.detect(image)
            if faces is None or len(faces) == 0:
                print(f"photo_positions: nenhum rosto em {name}; usará a posição padrão.")
                continue
            face = max(faces, key=lambda f: f[2] * f[3])
            position = position_for(width, height, face)
            if position is None:
                print(f"photo_positions: {name} é quadrada; sem ajuste.")
                continue
            results[name] = f"{position[0]}% {position[1]}%"
            print(f"photo_positions: {name} -> {results[name]}")
        except Exception as error:  # uma foto problemática não pode atrapalhar as outras
            print(f"photo_positions: erro em {name}: {str(error).splitlines()[0]}; mantido o valor anterior, se houver.")
            if name in previous:
                results[name] = previous[name]

    lines = ["# Gerado automaticamente por scripts/photo_positions.py a cada deploy. Não edite.",
             "# Para forçar um enquadramento, use photo_position em _data/people.yml."]
    if results:
        lines += [f'"{name}": "{value}"' for name, value in sorted(results.items())]
    else:
        lines.append("{}")
    tmp = OUT + ".tmp"
    with open(tmp, "w", encoding="utf8") as handle:
        handle.write("\n".join(lines) + "\n")
    os.replace(tmp, OUT)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:  # última rede de segurança
        print(f"photo_positions: falhou ({error}); nada foi alterado.")
        sys.exit(0)
