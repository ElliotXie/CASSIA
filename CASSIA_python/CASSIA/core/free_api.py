"""
Free API Access Module for CASSIA

Provides automatic free API key fetching for users without configured keys.
Keys are obtained from cassia.bio server with rate limiting.

Limits:
- 2 cluster annotations per machine (lifetime)

Supported providers:
- google (Google AI Studio / Gemini)
- together (Together AI)
- openrouter (OpenRouter)
"""

import os
import hashlib
import platform
import uuid
import socket
import threading
from typing import Optional, Dict, Tuple
from datetime import datetime

try:
    import requests
except ImportError:
    requests = None

# =============================================================================
# CONFIGURATION
# =============================================================================

# Server configuration
FREE_API_SERVER = os.environ.get("CASSIA_FREE_API_SERVER", "https://cassia.bio")
FREE_API_ENDPOINT = f"{FREE_API_SERVER}/api/free-api/get-key"
FREE_API_HEALTH_ENDPOINT = f"{FREE_API_SERVER}/api/free-api/health"
FREE_API_USAGE_ENDPOINT = f"{FREE_API_SERVER}/api/free-api/usage"

# Lifetime limit: 2 free cluster annotations per machine, total
MAX_FREE_CLUSTERS = 2

# Deprecated aliases (kept for backward compatibility)
MAX_CLUSTERS_PER_JOB = MAX_FREE_CLUSTERS
MAX_JOBS_PER_DAY = 1  # No longer used, kept to avoid import errors

# Supported providers for free API
FREE_API_PROVIDERS = {"google", "together", "openrouter"}

# =============================================================================
# THREAD-LOCAL STORAGE
# =============================================================================

# Thread-local storage for passing cluster count context between functions
_thread_local = threading.local()

# =============================================================================
# CACHING
# =============================================================================

# Thread-safe cache for session key
_free_key_cache: Dict[str, Dict] = {}
_cache_lock = threading.Lock()


def _get_cached_key(provider: str) -> Optional[str]:
    """Get cached free API key if still valid."""
    with _cache_lock:
        cache_key = f"free_api:{provider}"
        cached = _free_key_cache.get(cache_key)

        if cached:
            # Check if key is still valid (cached for current session)
            # Keys are cached until process ends or cache is cleared
            return cached.get('api_key')

        return None


def _cache_key(provider: str, api_key: str, remaining_clusters: int) -> None:
    """Cache a free API key for the session."""
    with _cache_lock:
        cache_key = f"free_api:{provider}"
        _free_key_cache[cache_key] = {
            'api_key': api_key,
            'remaining_clusters': remaining_clusters,
            'cached_at': datetime.utcnow().isoformat()
        }


def clear_free_api_cache() -> None:
    """Clear cached free API keys. Call this to force re-fetching keys."""
    with _cache_lock:
        _free_key_cache.clear()


# =============================================================================
# MACHINE FINGERPRINTING
# =============================================================================

def _generate_machine_id() -> str:
    """
    Generate a stable machine fingerprint for user identification.

    Uses multiple hardware/software identifiers to create a stable hash:
    - Platform info (OS, version)
    - Machine name (hostname)
    - MAC address (if available)
    - CPU info

    Returns:
        SHA-256 hash of combined identifiers (first 32 chars)

    Note:
        This is not meant to be cryptographically secure or impossible to spoof.
        It's a reasonable identifier for rate limiting free tier access.
    """
    components = []

    # Platform information
    try:
        components.append(platform.system())
        components.append(platform.release())
        components.append(platform.machine())
    except Exception:
        pass

    # Machine name (hostname)
    try:
        components.append(socket.gethostname())
    except Exception:
        pass

    # MAC address (uuid.getnode returns MAC-based value)
    try:
        mac = uuid.getnode()
        # Check if this is a real MAC (not random)
        # Random MACs have the multicast bit set
        if (mac >> 40) % 2 == 0:  # Not a random MAC
            components.append(str(mac))
    except Exception:
        pass

    # CPU/processor info
    try:
        proc = platform.processor()
        if proc:
            components.append(proc)
    except Exception:
        pass

    # Combine and hash
    combined = "|".join(str(c) for c in components if c)

    # Use SHA-256 and take first 32 characters
    machine_hash = hashlib.sha256(combined.encode('utf-8')).hexdigest()[:32]

    return machine_hash


# =============================================================================
# MAIN API FUNCTIONS
# =============================================================================

def get_remaining_free_clusters() -> int:
    """
    Query the server for how many free cluster annotations remain for this machine.

    Returns:
        Number of remaining free clusters (0 to MAX_FREE_CLUSTERS).
        Returns MAX_FREE_CLUSTERS if server is unreachable (optimistic default).
    """
    if requests is None:
        return MAX_FREE_CLUSTERS  # Optimistic: let server decide

    machine_id = _generate_machine_id()

    try:
        response = requests.get(
            FREE_API_USAGE_ENDPOINT,
            params={"machine_id": machine_id},
            timeout=5
        )
        if response.status_code == 200:
            data = response.json()
            return data.get('clusters_remaining', 0)
    except Exception:
        pass

    return MAX_FREE_CLUSTERS  # Optimistic default


def get_free_api_key(
    provider: str,
    num_clusters: int = 1
) -> Tuple[Optional[str], Optional[str]]:
    """
    Fetch a free API key from the CASSIA server.

    This function is called automatically by call_llm() when no API key is
    configured. It contacts the CASSIA server to obtain a temporary API key
    for the specified provider.

    Args:
        provider: The LLM provider ('google', 'together', 'openrouter')
        num_clusters: Number of clusters in this batch job (for rate limiting)

    Returns:
        Tuple of (api_key, error_message):
        - On success: (api_key_string, None)
        - On failure: (None, error_message_string)

    Rate Limits:
        - Maximum 2 cluster annotations per machine (lifetime)

    Example:
        >>> key, error = get_free_api_key("openrouter", num_clusters=2)
        >>> if error:
        ...     print(f"Failed: {error}")
        ... else:
        ...     # Use key for API call
    """
    # Check if requests is available
    if requests is None:
        return None, (
            "The 'requests' package is required for free API access. "
            "Install it with: pip install requests"
        )

    # Validate provider
    provider_lower = provider.lower()
    if provider_lower not in FREE_API_PROVIDERS:
        return None, (
            f"Free API access is not available for provider '{provider}'. "
            f"Supported providers: {', '.join(sorted(FREE_API_PROVIDERS))}. "
            f"Please set your own API key with CASSIA.set_api_key('{provider}', 'your-key')"
        )

    # Check cache first (avoid unnecessary server calls)
    cached_key = _get_cached_key(provider_lower)
    if cached_key:
        return cached_key, None

    # Generate machine ID
    machine_id = _generate_machine_id()

    # Fetch from server
    try:
        response = requests.post(
            FREE_API_ENDPOINT,
            json={
                "machine_id": machine_id,
                "provider": provider_lower,
                "num_clusters": num_clusters
            },
            headers={
                "Content-Type": "application/json",
                "User-Agent": "CASSIA-Python/1.0"
            },
            timeout=15
        )

        if response.status_code == 200:
            data = response.json()
            api_key = data.get('api_key')
            remaining_clusters = data.get('remaining_clusters', 0)

            if api_key:
                # Cache the key for this session
                _cache_key(provider_lower, api_key, remaining_clusters)

                # Print info message (only once per provider)
                print(f"[CASSIA] Using free API access for {provider_lower} "
                      f"({remaining_clusters} free cluster annotation(s) remaining)")

                return api_key, None
            else:
                return None, "Server returned invalid response (no api_key)"

        elif response.status_code == 429:
            # Rate limited (lifetime limit reached)
            try:
                data = response.json()
                error_msg = data.get('error', 'Free API limit reached')
                return None, (
                    f"{error_msg}\n"
                    f"To continue using CASSIA, set your own API key:\n"
                    f"  CASSIA.set_api_key('{provider}', 'your-key')"
                )
            except Exception:
                return None, (
                    f"Free API limit reached ({MAX_FREE_CLUSTERS} cluster annotations lifetime). "
                    f"Set your own API key: CASSIA.set_api_key('{provider}', 'your-key')"
                )

        elif response.status_code == 400:
            try:
                data = response.json()
                return None, data.get('error', 'Invalid request')
            except Exception:
                return None, 'Invalid request'

        elif response.status_code == 503:
            try:
                data = response.json()
                return None, data.get('error', 'Service temporarily unavailable')
            except Exception:
                return None, 'Free API service temporarily unavailable'

        else:
            return None, f"Server error: HTTP {response.status_code}"

    except requests.exceptions.Timeout:
        return None, (
            "Connection to CASSIA server timed out. "
            "Please check your internet connection and try again, "
            "or set your own API key."
        )
    except requests.exceptions.ConnectionError:
        return None, (
            "Could not connect to CASSIA server. "
            "Please check your internet connection and try again, "
            "or set your own API key."
        )
    except Exception as e:
        return None, f"Failed to fetch free API key: {str(e)}"


def is_free_api_available() -> bool:
    """
    Check if the free API server is reachable.

    Returns:
        True if server is reachable and responding, False otherwise
    """
    if requests is None:
        return False

    try:
        response = requests.get(FREE_API_HEALTH_ENDPOINT, timeout=5)
        if response.status_code == 200:
            data = response.json()
            return data.get('status') == 'ok'
        return False
    except Exception:
        return False


def get_free_api_usage() -> Optional[Dict]:
    """
    Get current lifetime usage statistics for this machine.

    Returns:
        Dictionary with usage info, or None if unavailable.

        Example return value:
        {
            'machine_id': 'abc12345...',
            'clusters_used': 1,
            'clusters_remaining': 1,
            'last_used_at': '2025-01-15T...',
            'limits': {
                'max_free_clusters': 2
            }
        }
    """
    if requests is None:
        return None

    machine_id = _generate_machine_id()

    try:
        response = requests.get(
            FREE_API_USAGE_ENDPOINT,
            params={"machine_id": machine_id},
            timeout=5
        )
        if response.status_code == 200:
            return response.json()
    except Exception:
        pass

    return None


def is_using_free_api(provider: str) -> bool:
    """
    Check if the given provider would use free API (no key set).

    Args:
        provider: Provider name to check

    Returns:
        True if no API key is set and provider supports free API
    """
    provider_lower = provider.lower()

    # Check if provider supports free API
    if provider_lower not in FREE_API_PROVIDERS:
        return False

    # Check if API key is set in environment
    env_var_map = {
        "openrouter": "OPENROUTER_API_KEY",
        "google": "GOOGLE_API_KEY",
        "together": "TOGETHER_API_KEY",
    }

    env_var = env_var_map.get(provider_lower)
    if env_var and os.environ.get(env_var):
        return False  # Key is set, not using free API

    return True


def get_free_api_info() -> Dict:
    """
    Get information about free API access.

    Returns:
        Dictionary with free API configuration and limits
    """
    return {
        'server': FREE_API_SERVER,
        'supported_providers': list(FREE_API_PROVIDERS),
        'limits': {
            'max_free_clusters': MAX_FREE_CLUSTERS
        },
        'is_available': is_free_api_available()
    }
