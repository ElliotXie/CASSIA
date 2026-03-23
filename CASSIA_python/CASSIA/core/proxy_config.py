"""
Proxy configuration for CASSIA API calls.

Allows users in regions where US API providers are blocked (e.g., China)
to route traffic through a Cloudflare Worker reverse proxy.
Users still provide their own API keys (BYOK) — only network traffic is relayed.

Auto-detection: On first import, CASSIA detects if the user is in China
(via timezone) and enables the proxy automatically. Users can override with:
    CASSIA.set_proxy(None)      # Disable proxy
    CASSIA.set_proxy("china")   # Re-enable proxy
"""

import time
import os

# Proxy endpoint mappings: proxy_name -> {provider -> proxied_base_url}
PROXY_ENDPOINTS = {
    "china": {
        "openai": "https://cassia-proxy.solodev66.workers.dev/openai/v1",
        "anthropic": "https://cassia-proxy.solodev66.workers.dev/anthropic",
        "openrouter": "https://cassia-proxy.solodev66.workers.dev/openrouter/api/v1",
    }
}

# Module-level active proxy state
# None = not yet initialized, False = explicitly disabled, str = active proxy name
_active_proxy = None
_auto_detected = False


def _detect_china():
    """
    Detect if the user is likely in China based on system timezone.
    Uses UTC offset (UTC+8) and timezone name heuristics.
    No network calls — purely local check.
    """
    try:
        # Check timezone name first (most reliable)
        tz_name = time.tzname[0] if time.tzname else ""
        china_tz_names = {"CST", "China Standard Time", "Asia/Shanghai", "Asia/Chongqing"}
        if tz_name in china_tz_names:
            return True

        # Check UTC offset: China is UTC+8 (28800 seconds)
        # time.timezone is seconds WEST of UTC (negative for east), and doesn't account for DST
        # UTC+8 means time.timezone == -28800
        utc_offset = -time.timezone
        if utc_offset == 28800:
            # UTC+8 could also be Singapore, Philippines, etc.
            # Check locale for additional signal
            import locale
            loc = locale.getdefaultlocale()[0] or ""
            if loc.startswith("zh") or loc.startswith("cn"):
                return True

        # Check LANG environment variable
        lang = os.environ.get("LANG", "") + os.environ.get("LC_ALL", "")
        if "zh_CN" in lang or "zh_TW" in lang:
            return True

    except Exception:
        pass

    return False


def _auto_detect_proxy():
    """Auto-detect and set proxy on first use. Called lazily."""
    global _active_proxy, _auto_detected
    if _auto_detected:
        return
    _auto_detected = True

    # Don't auto-detect if user has explicitly set CASSIA_PROXY env var
    env_proxy = os.environ.get("CASSIA_PROXY")
    if env_proxy is not None:
        if env_proxy.lower() in ("", "none", "off", "false", "0"):
            _active_proxy = False
        else:
            _active_proxy = env_proxy.lower().strip()
        return

    if _detect_china():
        _active_proxy = "china"


def set_proxy(proxy_name):
    """
    Set the active proxy for all CASSIA API calls.

    Args:
        proxy_name (str or None or False):
            - "china": Route through China proxy
            - None or False: Disable proxy (direct connection)

    Raises:
        ValueError: If proxy_name is not a recognized preset.
    """
    global _active_proxy, _auto_detected
    _auto_detected = True  # User has explicitly set, skip auto-detection

    if proxy_name is None or proxy_name is False:
        _active_proxy = False
        return
    proxy_name = str(proxy_name).lower().strip()
    if proxy_name in ("none", "off", "false"):
        _active_proxy = False
        return
    if proxy_name not in PROXY_ENDPOINTS:
        raise ValueError(
            f"Unknown proxy '{proxy_name}'. "
            f"Available proxies: {', '.join(sorted(PROXY_ENDPOINTS.keys()))}"
        )
    _active_proxy = proxy_name


def get_proxy():
    """Return the currently active proxy name, or None if disabled."""
    _auto_detect_proxy()
    if _active_proxy is False:
        return None
    return _active_proxy


def resolve_proxy_url(provider):
    """
    Resolve the proxied base URL for a given provider.

    Args:
        provider (str): The provider name (e.g., "openai", "anthropic", "openrouter").

    Returns:
        str or None: The proxied base URL if a proxy is active and the provider is supported,
                     otherwise None (meaning use the default endpoint).
    """
    _auto_detect_proxy()
    if not _active_proxy or _active_proxy is False:
        return None
    endpoints = PROXY_ENDPOINTS.get(_active_proxy, {})
    return endpoints.get(provider.lower())
