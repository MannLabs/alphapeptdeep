import torch


def _is_mps_available() -> bool:
    try:
        return torch.backends.mps.is_available()
    except AttributeError:
        return False


torch_devices: dict = {
    "gpu": {
        "is_available": torch.cuda.is_available,
        "device": "cuda",
    },
    "cuda": {
        "is_available": torch.cuda.is_available,
        "device": "cuda",
    },
    "mps": {
        "is_available": _is_mps_available,
        "device": "mps",
    },
    "m1": {
        "is_available": _is_mps_available,
        "device": "mps",
    },
}


def get_device(device: str, device_ids: list = []) -> tuple:
    """Device name to torch.device

    Parameters
    ----------
    device : str
        device type: cuda, gpu, mps, m1.
    device_ids : list, optional
        device_ids for cuda, by default []

    Returns
    -------
    tuple
        torch.device
        str: device name
    """
    device = device.lower()
    if device in torch_devices:
        if torch_devices[device]["is_available"]():
            if torch_devices[device]["device"] == "cuda" and len(device_ids) > 0:
                return torch.device(
                    f'cuda:{",".join(str(_id) for _id in device_ids)}'
                ), "cuda"
            else:
                return (torch.device(torch_devices[device]["device"]), device)
    return torch.device("cpu"), "cpu"


def get_available_device() -> tuple:
    for name, item in torch_devices.items():
        if item["is_available"]():
            return torch.device(item["device"]), name
    return torch.device("cpu"), "cpu"
