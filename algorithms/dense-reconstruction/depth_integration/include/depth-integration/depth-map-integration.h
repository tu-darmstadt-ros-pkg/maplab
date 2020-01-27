#ifndef DEPTH_INTEGRATION_DEPTH_MAP_INTEGRATION_H_
#define DEPTH_INTEGRATION_DEPTH_MAP_INTEGRATION_H_

#include <functional>
#include <unordered_set>

#include <vi-map/vi-map.h>
#include <voxblox/core/common.h>

namespace depth_integration {

typedef std::function<void(
    const aslam::Transformation& /*T_G_C*/, const cv::Mat& /*depth_map*/,
    const cv::Mat& /*intensities*/, const aslam::Camera& /*camera*/)>
    DepthMapIntegrationFunction;

// This integration function will provide one transformation for every line of
// the depth map, based on the line-delay defined in the camera. This allows us
// to either integrate every line using a different transformation (motion
// undistortion) or just to simply create an undistorted depth map.
typedef std::function<void(
    const aslam::TransformationVector& /*T_G_C_vec*/,
    const cv::Mat& /*depth_map*/, const cv::Mat& /*intensities*/,
    const aslam::Camera& /*camera*/)>
    DepthMapUndistortionAndIntegrationFunction;

static std::unordered_set<backend::ResourceType, backend::ResourceTypeHash>
    kSupportedDepthMapInputTypes{backend::ResourceType::kRawDepthMap,
                                 backend::ResourceType::kOptimizedDepthMap};

// Calls the integration function for all depth (optional and frame) resources.
template <typename DepthIntegrationFuctionType>
void integrateAllDepthMapResourcesOfType(
    const vi_map::MissionIdList& mission_ids,
    const backend::ResourceType& input_resource_type,
    const vi_map::VIMap& vi_map,
    DepthIntegrationFuctionType integration_function);

// Calls the integration function for all depth frame resources from the
// selected missions using the integration function
template <typename DepthIntegrationFuctionType>
void integrateAllFrameDepthMapResourcesOfType(
    const vi_map::MissionIdList& mission_ids,
    const backend::ResourceType& input_resource_type,
    const vi_map::VIMap& vi_map,
    DepthIntegrationFuctionType integration_function);

// Calls the integration function for all optional depth resources from the
// selected missions using the integration function.
template <typename DepthIntegrationFuctionType>
void integrateAllOptionalSensorDepthMapResourcesOfType(
    const vi_map::MissionIdList& mission_ids,
    const backend::ResourceType& input_resource_type,
    const vi_map::VIMap& vi_map,
    DepthIntegrationFuctionType integration_function);

// Calls the integration function for all depth (optional and frame) resources.
template <typename DepthIntegrationFuctionType>
void integrateAllDepthMapResourcesOfType(
    const vi_map::MissionIdList& mission_ids,
    const backend::ResourceType& input_resource_type,
    const vi_map::VIMap& vi_map,
    DepthIntegrationFuctionType integration_function) {
  CHECK(integration_function);

  // Integrate all frame resources.
  integrateAllFrameDepthMapResourcesOfType(
      mission_ids, input_resource_type, vi_map, integration_function);

  // Integrate all optional sensor resources.
  integrateAllOptionalSensorDepthMapResourcesOfType(
      mission_ids, input_resource_type, vi_map, integration_function);
}

}  // namespace depth_integration

#endif  // DEPTH_INTEGRATION_DEPTH_MAP_INTEGRATION_H_
