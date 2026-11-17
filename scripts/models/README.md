# Face detection model

`face_detection_yunet_2023mar.onnx` is the **YuNet** face detector from the
[OpenCV Zoo](https://github.com/opencv/opencv_zoo/tree/main/models/face_detection_yunet) (MIT license).
It is used by `scripts/photo_positions.py` to find where the face is in each photo of `assets/img/people/`.
The model only returns the position of the face (a box); no photo leaves the machine that runs the script.
