"""
Proxy configuration for CASSIA API calls.

Allows users in regions where US API providers are blocked (e.g., China)
to route traffic through a Cloudflare Worker reverse proxy.
Users still provide their own API keys (BYOK) — only network traffic is relayed.

EMERGENCY KILL-SWITCH (2026-04-07):
Auto-detection is currently DISABLED. The deployed Cloudflare worker lives at
`cassia-proxy.solodev66.workers.dev` and the `*.workers.dev` apex is DNS-poisoned
in mainland China, so routing China users through it just makes Anthropic /
OpenRouter (which would otherwise work directly) fail too. The previous
`time.tzname == "CST"` check also collided with US Central Standard Time and
mis-routed every Chicago/Madison/Texas user through the same broken proxy.

The proxy now only activates when the user *explicitly* opts in via:
    CASSIA.set_proxy("china")              # in code
    export CASSIA_PROXY=china               # via environment

To re-enable auto-detection later (after the worker is rebound to a custom
domain), restore the `_detect_china()` call in `_auto_detect_proxy()`.
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

    NOTE: Currently NOT called from `_auto_detect_proxy()` (see kill-switch
    note at the top of this file). Kept here so it can be re-enabled later
    once the upstream worker is moved off `*.workers.dev`.

    The previous version matched on `time.tzname[0] == "CST"`, which is
    ambiguous: "CST" is also the abbreviation for US Central Standard Time
    (Chicago/Madison/Texas). That bug mis-classified every US-Central user
    as a China user. This version requires an unambiguous IANA-style name
    or an explicit `zh_CN` / `zh_TW` LANG variable.
    """
    try:
        # Only trust unambiguous IANA-style China timezone names.
        # Notably absent: "CST" (collides with US Central Standard Time).
        tz_name = time.tzname[0] if time.tzname else ""
        china_tz_names = {
            "China Standard Time",
            "Asia/Shanghai",
            "Asia/Chongqing",
            "Asia/Harbin",
            "Asia/Urumqi",
            "PRC",
        }
        if tz_name in china_tz_names:
            return True

        # LANG / LC_ALL set explicitly to a mainland-China locale.
        # We deliberately do NOT use UTC+8 + zh_* — that mis-classifies
        # users in Singapore, Hong Kong, Taiwan, and overseas Chinese.
        lang = os.environ.get("LANG", "") + os.environ.get("LC_ALL", "")
        if "zh_CN" in lang:
            return True

    except Exception:
        pass

    return False


def _auto_detect_proxy():
    """
    Resolve the proxy from explicit user opt-in only.

    Auto-detection of mainland China users is INTENTIONALLY DISABLED — see
    the kill-switch note at the top of this file. This function now only
    activates a proxy when the user has explicitly opted in via the
    `CASSIA_PROXY` environment variable or a `set_proxy()` call.
    """
    global _active_proxy, _auto_detected
    if _auto_detected:
        return
    _auto_detected = True

    # Honor explicit user opt-in via environment variable.
    env_proxy = os.environ.get("CASSIA_PROXY")
    if env_proxy is not None:
        if env_proxy.lower() in ("", "none", "off", "false", "0"):
            _active_proxy = False
        else:
            _active_proxy = env_proxy.lower().strip()
        return

    # Auto-detection intentionally disabled — see kill-switch note above.
    # (Previously: `if _detect_china(): _active_proxy = "china"`)
    _active_proxy = False


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
